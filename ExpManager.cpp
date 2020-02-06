// ***************************************************************************************************************
//
//          Mini-Aevol is a reduced version of Aevol -- An in silico experimental evolution platform
//
// ***************************************************************************************************************
//
// Copyright: See the AUTHORS file provided with the package or <https://gitlab.inria.fr/rouzaudc/mini-aevol>
// Web: https://gitlab.inria.fr/rouzaudc/mini-aevol
// E-mail: See <jonathan.rouzaud-cornabas@inria.fr>
// Original Authors : Jonathan Rouzaud-Cornabas
//
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ***************************************************************************************************************


#include <cmath>
#include <map>
#include <sys/stat.h>
#include <err.h>
#include <chrono>
#include <iostream>

#ifdef USE_CUDA
#include "nvToolsExt.h"
#include <cuda_profiler_api.h>
using namespace std::chrono;
#endif

using namespace std;

#include "ExpManager.h"
#include "Algorithms.h"
#include "AeTime.h"
#include "RNA.h"
#include "Protein.h"
#include "Organism.h"
#include "Gaussian.h"

#include <utility>

/**
 * Constructor for initializing a new simulation
 *
 * @param grid_height : Height of the grid containing the organisms
 * @param grid_width : Width of the grid containing the organisms
 * @param seed : Global seed for all the PRNG of the simulation
 * @param mutation_rate : Mutation rates for all the organism during the simulation
 * @param init_length_dna : Size of the randomly generated DNA at the initialization of the simulation
 * @param w_max : Maximum width of the triangle generated by a Protein
 * @param selection_pressure : Selection pressure used during the selection process
 * @param backup_step : How much often checkpoint must be done
 */
ExpManager::ExpManager(int grid_height, int grid_width, int seed, double mutation_rate, int init_length_dna,
                       double w_max, int selection_pressure, int backup_step)
: rng_(grid_width, grid_height, seed)
{
    // Initializing the data structure
    grid_height_ = grid_height;
    grid_width_ = grid_width;

    backup_step_ = backup_step;

    nb_indivs_ = grid_height * grid_width;

    w_max_ = w_max;
    selection_pressure_ = selection_pressure;

    internal_organisms_ = new Organism*[nb_indivs_];
    prev_internal_organisms_ = new Organism*[nb_indivs_];
    internal_ids_ = new OrganismIDs[nb_indivs_];
    prev_internal_ids_ = new OrganismIDs[nb_indivs_];

    next_generation_reproducer_ = new int[nb_indivs_]();

    mutation_rate_ = mutation_rate;

    // Building the target environment
    Gaussian *g1 = new Gaussian(1.2, 0.52, 0.12);
    Gaussian *g2 = new Gaussian(-1.4, 0.5, 0.07);
    Gaussian *g3 = new Gaussian(0.3, 0.8, 0.03);

    target = new double[300];
    for (int i = 0; i < 300; i++) {
        double pt_i = ((double) i) / 300.0;

        double tmp = g1->compute_y(pt_i);
        tmp += g2->compute_y(pt_i);
        tmp += g3->compute_y(pt_i);

        tmp = tmp > Y_MAX ? Y_MAX : tmp;
        tmp = tmp < Y_MIN ? Y_MIN : tmp;

        target[i] = tmp;
    }

    delete g1;
    delete g2;
    delete g3;

    geometric_area_ = 0;


    for (int i = 0; i < 299; i++) {
        geometric_area_ += ((fabs(target[i]) + fabs(target[i + 1])) / (600.0));
    }

    printf("Initialized environmental target %f\n", geometric_area_);

    // Generate a random organism that is better than nothing
    Organism* parent;
    for (;;) {
        Organism indiv(rng_.gen(0, Threefry::MUTATION), init_length_dna);

        indiv.start_stop_RNA();
        indiv.compute_RNA();

        //indiv.start_protein();
        //indiv.compute_protein();
        indiv.compute_proteins();

        indiv.translate_protein(w_max);

        indiv.compute_phenotype_fitness(selection_pressure, target);

        double r_compare = round((indiv.metaerror - geometric_area_) * 1E10) / 1E10;
        if (!(r_compare >= 0))
        {
            parent = new Organism(std::move(indiv), 0);
            break;
        }
    }

    printf("Populating the environment\n");
    int parent_length = parent->dna_.seq_.size();

    int global_id = AeTime::time() * nb_indivs_;
    // Create a population of clones based on the randomly generated organism

    #pragma omp parallel for
    for (int indiv_id = 0; indiv_id < nb_indivs_; indiv_id++) {
        prev_internal_organisms_[indiv_id] = internal_organisms_[indiv_id] = new Organism(*parent, 0);

        OrganismIDs& ids = internal_ids_[indiv_id];
        ids.indiv_id_ = indiv_id;
        ids.parent_id_ = 0;
        ids.global_id_ = global_id++;
        ids.parent_length_ = parent_length;
        prev_internal_ids_[indiv_id] = ids;
    }

    delete parent;

    // Create backup and stats directory
    create_directory();
}

/**
 * Constructor to resume/restore a simulation from a given backup/checkpoint file
 *
 * @param time : resume from this generation
 */
ExpManager::ExpManager(int time)
{
    target = new double[300];

    load(time);

    geometric_area_ = 0;
    for (int i = 0; i < 299; i++) {
        geometric_area_ += ((fabs(target[i]) + fabs(target[i + 1])) / (600.0));
    }

    printf("Initialized environmental target %f\n", geometric_area_);
}

/**
 * Create stats and backup directory
 */
void ExpManager::create_directory() {
    // Backup
    int status = mkdir("backup", 0755);
    if (status == -1 && errno != EEXIST) {
        err(EXIT_FAILURE, "backup");
    }

    // Stats
    status = mkdir("stats", 0755);
    if (status == -1 && errno != EEXIST) {
        err(EXIT_FAILURE, "stats");
    }
}

/**
 * Checkpointing/Backup of the population of organisms
 *
 * @param t : simulated time of the checkpoint
 */
void ExpManager::save(int t) const {
    char exp_backup_file_name[255];

    sprintf(exp_backup_file_name, "backup/backup_%d.zae", t);

    // -------------------------------------------------------------------------
    // Open backup files
    // -------------------------------------------------------------------------
    gzFile exp_backup_file = gzopen(exp_backup_file_name, "w");


    // -------------------------------------------------------------------------
    // Check that files were correctly opened
    // -------------------------------------------------------------------------
    if (exp_backup_file == Z_NULL) {
        printf("Error: could not open backup file %s\n",
               exp_backup_file_name);
        exit(EXIT_FAILURE);
    }


    // -------------------------------------------------------------------------
    // Write the backup file
    // -------------------------------------------------------------------------
    gzwrite(exp_backup_file, &t, sizeof(t));

    gzwrite(exp_backup_file, &grid_height_, sizeof(grid_height_));
    gzwrite(exp_backup_file, &grid_width_, sizeof(grid_width_));

    gzwrite(exp_backup_file, &nb_indivs_, sizeof(nb_indivs_));

    gzwrite(exp_backup_file, &backup_step_, sizeof(backup_step_));

    gzwrite(exp_backup_file, &w_max_, sizeof(w_max_));
    gzwrite(exp_backup_file, &selection_pressure_, sizeof(selection_pressure_));

    gzwrite(exp_backup_file, &mutation_rate_, sizeof(mutation_rate_));

    /*
    for (int i = 0; i < 300; i++) {
        double tmp = target[i];
        gzwrite(exp_backup_file, &tmp, sizeof(tmp));
    }
    */
    gzwrite(exp_backup_file, &target[0], 300 * sizeof(target[0]));

    for (int indiv_id = 0; indiv_id < nb_indivs_; indiv_id++) {
        OrganismIDs& ids = prev_internal_ids_[indiv_id];
        ids.indiv_id_ = ids.global_id_ % nb_indivs_;
        gzwrite(exp_backup_file, &ids, sizeof(ids));
        prev_internal_organisms_[indiv_id]->save(exp_backup_file);
    }

    rng_.save(exp_backup_file);

    if (gzclose(exp_backup_file) != Z_OK) {
        cerr << "Error while closing backup file" << endl;
    }
}

/**
 * Loading a simulation from a checkpoint/backup file
 *
 * @param t : resuming the simulation at this generation
 */
void ExpManager::load(int t) {
    char exp_backup_file_name[255];

    sprintf(exp_backup_file_name, "backup/backup_%d.zae", t);

    // -------------------------------------------------------------------------
    // Open backup files
    // -------------------------------------------------------------------------
    gzFile exp_backup_file = gzopen(exp_backup_file_name, "r");


    // -------------------------------------------------------------------------
    // Check that files were correctly opened
    // -------------------------------------------------------------------------
    if (exp_backup_file == Z_NULL) {
        printf("Error: could not open backup file %s\n",
               exp_backup_file_name);
        exit(EXIT_FAILURE);
    }


    // -------------------------------------------------------------------------
    // Write the backup file
    // -------------------------------------------------------------------------
    int time;
    gzread(exp_backup_file, &time, sizeof(time));
    AeTime::set_time(time);

    gzread(exp_backup_file, &grid_height_, sizeof(grid_height_));

    gzread(exp_backup_file, &grid_width_, sizeof(grid_width_));

    gzread(exp_backup_file, &nb_indivs_, sizeof(nb_indivs_));

    internal_organisms_ = new Organism*[nb_indivs_];
    prev_internal_organisms_ = new Organism*[nb_indivs_];
    internal_ids_ = new OrganismIDs[nb_indivs_];
    prev_internal_ids_ = new OrganismIDs[nb_indivs_];

    // No need to save/load this field from the backup because it will be set at selection()
    next_generation_reproducer_ = new int[nb_indivs_]();

    gzread(exp_backup_file, &backup_step_, sizeof(backup_step_));

    gzread(exp_backup_file, &w_max_, sizeof(w_max_));
    gzread(exp_backup_file, &selection_pressure_, sizeof(selection_pressure_));

    gzread(exp_backup_file, &mutation_rate_, sizeof(mutation_rate_));

    for (int i = 0; i < 300; i++) {
        double tmp;
        gzread(exp_backup_file, &tmp, sizeof(tmp));
        target[i] = tmp;
    }

    for (int indiv_id = 0; indiv_id < nb_indivs_; indiv_id++) {
        OrganismIDs& ids = internal_ids_[indiv_id];
        gzread(exp_backup_file, &ids, sizeof(ids));
        prev_internal_ids_[indiv_id] = internal_ids_[indiv_id];
        prev_internal_organisms_[indiv_id] = internal_organisms_[indiv_id] = new Organism(exp_backup_file);
        // promoters have to be recomputed, they are not save in the backup
        internal_organisms_[indiv_id]->start_stop_RNA();
    }

    rng_ = Threefry(grid_width_, grid_height_, exp_backup_file);

    if (gzclose(exp_backup_file) != Z_OK) {
        cerr << "Error while closing backup file" << endl;
    }
}

/**
 * Prepare the mutation generation of an organism
 *
 * @param indiv_id : Organism unique id
 */
bool ExpManager::prepare_mutation(int indiv_id) {
    int parent_id = next_generation_reproducer_[indiv_id];
    Organism* parent = prev_internal_organisms_[parent_id];
    int parent_length = parent->dna_.seq_.size();

    DnaMutator dna_mutator(
            Threefry::Gen(rng_.gen(indiv_id, Threefry::MUTATION)),
            parent_length, mutation_rate_
    );

    bool mutations = dna_mutator.hasMutate();
    if (mutations) {
        Organism* child = new Organism(*parent, 0);
        internal_organisms_[indiv_id] = child;

        OrganismIDs& ids = internal_ids_[indiv_id];
        ids.indiv_id_ = prev_internal_ids_[parent_id].indiv_id_;
        ids.parent_id_ = parent_id;
        ids.global_id_ = AeTime::time() * nb_indivs_ + indiv_id;
        ids.parent_length_ = parent_length;

        dna_mutator.apply_mutations(*child);
    } else {
        internal_organisms_[indiv_id] = parent;
        internal_ids_[indiv_id] = prev_internal_ids_[parent_id];

        parent->usage_count_++;
        parent->reset_mutation_stats();
    }

    return mutations;
}

/**
 * Destructor of the ExpManager class
 */
ExpManager::~ExpManager() {
    delete[] internal_organisms_;

    for (auto i = 0; i < nb_indivs_; ++i) {
        if (--prev_internal_organisms_[i]->usage_count_ == 0) delete prev_internal_organisms_[i];
    }
    delete[] prev_internal_organisms_;
    delete[] internal_ids_;
    delete[] prev_internal_ids_;
    delete[] next_generation_reproducer_;
    delete[] target;
}

/**
 * Execute a generation of the simulation for all the Organisms
 *
 * @param w_max : Maximum width of the triangle generated by a Protein
 * @param selection_pressure : Selection pressure used during the selection process
 * @param first_gen : is it the first generation simulated ? (generation 1 or first generation after a restore)
 */
void ExpManager::run_a_step(double w_max, double selection_pressure, bool first_gen) {

    // Running the simulation process for each organism
    double best_fitness;
    int idx_best = 0;

    #pragma omp parallel
    {
        #pragma omp for
        for (int indiv_id = 0; indiv_id < nb_indivs_; indiv_id++) {
            selection(indiv_id);
            if (prepare_mutation(indiv_id)) {
                Organism &indiv = (*internal_organisms_[indiv_id]);
                indiv.opt_prom_compute_RNA();
                indiv.compute_proteins();
                indiv.translate_protein(w_max);
                indiv.compute_phenotype_fitness(selection_pressure, target);
                indiv.compute_protein_stats();
            }
        }

        #pragma omp for
        for (int indiv_id = 0; indiv_id < nb_indivs_; indiv_id++) {
            if (--prev_internal_organisms_[indiv_id]->usage_count_ == 0) delete prev_internal_organisms_[indiv_id];

            prev_internal_organisms_[indiv_id] = internal_organisms_[indiv_id];
            internal_organisms_[indiv_id] = nullptr;
        }

        #pragma omp single
        {
            OrganismIDs *tmp = prev_internal_ids_;
            prev_internal_ids_ = internal_ids_;
            internal_ids_ = tmp;

            // Search for the best
            best_fitness = prev_internal_organisms_[0]->fitness;
        }

        #pragma omp for
        for (int indiv_id = 1; indiv_id < nb_indivs_; indiv_id++) {
            if (prev_internal_organisms_[indiv_id]->fitness > best_fitness) {
                idx_best = indiv_id;
                best_fitness = prev_internal_organisms_[indiv_id]->fitness;
            }
        }
    }

    best_indiv = prev_internal_organisms_[idx_best];
    // Stats

    stats_best.reinit(AeTime::time(), true);
    stats_mean.reinit(AeTime::time(), false);

    if (first_gen) {
        stats_best.prepare();
        stats_mean.prepare();

		for (int indiv_id = 0; indiv_id < nb_indivs_; indiv_id++) {
		    prev_internal_organisms_[indiv_id]->compute_protein_stats();
		}
    }

    stats_best.compute_best(*(best_indiv));
    stats_mean.compute_average(prev_internal_organisms_, nb_indivs_);

    stats_best.write_best();
    stats_mean.write_average();
}


/**
 * Selection process: for a given cell in the grid of the population, compute which organism win the computation
 *
  * @param indiv_id : Unique identification number of the cell
 */
void ExpManager::selection(int indiv_id) {
    int8_t selection_scope_x = 3;
    int8_t selection_scope_y = 3;
    int8_t neighborhood_size = 9;

    double local_fit_array[neighborhood_size];
    double probs[neighborhood_size];
    int count = 0;
    double sum_local_fit = 0.0;

    int32_t x = indiv_id / grid_height_;
    int32_t y = indiv_id % grid_height_;

    int cur_x, cur_y;

    for (int8_t i = -1; i < selection_scope_x - 1; i++) {
        #pragma omp simd
        for (int8_t j = -1; j < selection_scope_y - 1; j++) {
            cur_x = (x + i + grid_width_) % grid_width_;
            cur_y = (y + j + grid_height_) % grid_height_;

            local_fit_array[count] = prev_internal_organisms_[cur_x * grid_height_ + cur_y]->fitness;

            sum_local_fit += local_fit_array[count];

            count++;
        }
    }

    for (int8_t i = 0; i < neighborhood_size; i++) {
        probs[i] = local_fit_array[i] / sum_local_fit;
    }

    Threefry::Gen rng = rng_.gen(indiv_id, Threefry::REPROD);
    int found_org = rng.roulette_random(probs, neighborhood_size);

    int x_offset = (found_org / selection_scope_x) - 1;
    int y_offset = (found_org % selection_scope_y) - 1;

    next_generation_reproducer_[indiv_id] = ((x + x_offset + grid_width_) % grid_width_) * grid_height_ +
                                            ((y + y_offset + grid_height_) % grid_height_);
}

/**
 * Run the evolution for a given number of generation
 *
 * @param nb_gen : Number of generations to simulate
 */
void ExpManager::run_evolution(int nb_gen) {

    #pragma omp parallel for
    for (int indiv_id = 0; indiv_id < nb_indivs_; indiv_id++) {
        Organism& indiv = (*internal_organisms_[indiv_id]);
        
        indiv.opt_prom_compute_RNA();

        indiv.compute_proteins();

        indiv.translate_protein(w_max_);

        indiv.compute_phenotype_fitness(selection_pressure_, target);
    }

    printf("Running evolution from %d to %d\n", AeTime::time(), AeTime::time() + nb_gen);
    bool firstGen = true;

    for (int gen = 0; gen < nb_gen; gen++) {
        AeTime::plusplus();

        run_a_step(w_max_, selection_pressure_, firstGen);

        firstGen = false;
        //printf("Generation %d : Best individual fitness %e\n", AeTime::time(), best_indiv->fitness);

        if (AeTime::time() % backup_step_ == 0) {
            save(AeTime::time());
            cout << "Backup for generation " << AeTime::time() << " done !" << endl;
        }
    }
}

#ifdef USE_CUDA
void ExpManager::run_evolution_on_gpu(int nb_gen) {
  cudaProfilerStart();
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  cout << "Transfer" << endl;
  transfer_in(this, true);
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration_transfer_in = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
  cout << "Transfer done in " << duration_transfer_in << endl;

    for (int indiv_id = 0; indiv_id < nb_indivs_; indiv_id++) {
        Organism& indiv = (*internal_organisms_[indiv_id]);

        /*
        DnaMutator dna_mutator(
                rng_.gen(indiv_id, Threefry::MUTATION),
                prev_internal_organisms_[next_generation_reproducer_[indiv_id]]->dna_.seq_.size(),
                mutation_rate_, indiv_id);
        */

        indiv.opt_prom_compute_RNA();
        //indiv.compute_RNA();

        //indiv.start_protein();
        //indiv.compute_protein();
        indiv.compute_proteins();

        indiv.translate_protein(w_max_);

        indiv.compute_phenotype_fitness(selection_pressure_, target);
    }

  printf("Running evolution GPU from %d to %d\n",AeTime::time(),AeTime::time()+nb_gen);
  bool firstGen = true;
  for (int gen = 0; gen < nb_gen+1; gen++) {
    if(gen == 91) nvtxRangePushA("generation 91 to 100");
    AeTime::plusplus();
      //run_a_step(w_max_,selection_pressure_,firstGen);

      high_resolution_clock::time_point t1 = high_resolution_clock::now();
      run_a_step_on_GPU(nb_indivs_, w_max_, selection_pressure_, grid_width_, grid_height_,mutation_rate_);

      t2 = high_resolution_clock::now();
      auto duration_transfer_in = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

      std::cout<<"LOG,"<<duration_transfer_in<<std::endl;

    firstGen = false;
    if(gen == 100) nvtxRangePop();
    printf("Generation %d : \n",AeTime::time());
  }
  cudaProfilerStop();
}
#endif
