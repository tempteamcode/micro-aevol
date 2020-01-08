#include <vector>

template <typename item_t>
class uninitialized_vector
{
public:
  uninitialized_vector(size_t size = 0) : used(size), allocated(size) { if (allocated > 0) data = new item_t[allocated]; }
  ~uninitialized_vector() { if (allocated > 0) delete[] data; }
  
  uninitialized_vector (const uninitialized_vector& other) : used(other.size), allocated(other.size)
  {
    if (allocated > 0) {
      data = new item_t[allocated];
      
      if (std::is_trivially_copyable<item_t>::value) {
        std::memcpy(data, other.data, used);
      } else {
        item_t* dest = data;
        item_t* src = other.data;
        size_t count = used;
        while (count --> 0) *dest++ = *src++;
      }
    }
  }
  
  uninitialized_vector (uninitialized_vector&& other) : used(other.used), allocated(other.allocated), data(other.data)
  {
    other.used = 0;
    other.allocated = 0;
  }
  
  item_t& operator[](size_t pos) { return data[pos]; }
  const item_t& operator[](size_t pos) const { return data[pos]; }
  
  void redim_discard(size_t size)
  {
    if (size > allocated)
    {
      if (allocated > 0) delete[] data;
      allocated = size;
      data = new item_t[size];
    }
    used = size;
  }
  
  void redim_preserve(size_t size)
  {
    if (size > allocated)
    {
      item_t* data_old = data;
      data = new item_t[size];
      
      if (allocated > 0) {
        if (used > 0) {
          if (std::is_trivially_copyable<item_t>::value) {
            std::memcpy(data, data_old, used);
          } else {
            item_t* dest = data;
            item_t* src = data_old;
            size_t count = used;
            while (count --> 0) *dest++ = *src++;
          }
        }
        
        delete[] data_old;
      }
      
      allocated = size;
    }
    used = size;
  }
  
  size_t used;
  size_t allocated;
  item_t* data;
};

class better_bitset
{
public:
  typedef unsigned int int_t;
  typedef uninitialized_vector<int_t> vec_t;
  
  static constexpr size_t sizeof_int = sizeof(int_t);
  
  better_bitset() = default;
  ~better_bitset() = default;
  
  inline better_bitset(size_t size) : used(size / sizeof_int + 1), data(used) { }
  
  inline size_t size() const { return used; }
  
  template<typename function_t>
  void generate(size_t size, function_t generator)
  {
    used = (size / sizeof_int) + 1;
    data.redim_discard(used);
  }
  
  
  const bool operator[](size_t pos)
  {
    throw "better_bitset operator[]"; // TODO
  }
  
  bool test(size_t pos);
  void set(size_t pos);
  void set(size_t pos, bool value);
  void reset(size_t pos);
  
  
  void range_erase(size_t begin, size_t end)
  {
    throw "better_bitset range_erase"; // TODO
  }
  
  void range_move(size_t begin, size_t end, size_t pos)
  {
    throw "better_bitset range_move"; // TODO
  }
  
  void range_insert(const better_bitset&, size_t pos)
  {
    throw "better_bitset range_insert"; // TODO
  }
  
private:
  size_t used;
  vec_t data;
};

