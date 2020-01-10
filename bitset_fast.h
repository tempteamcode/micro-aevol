#include <vector>
#include <utility>

template <typename item_t>
class uninitialized_vector
{
public:
  uninitialized_vector(size_t size = 0) : used(size), allocated(size) { if (allocated > 0) data = new item_t[allocated]; }
  ~uninitialized_vector() { if (allocated > 0) delete[] data; }
  
  uninitialized_vector& operator=(const uninitialized_vector& other)
  {
    if (&other == this) return *this;
    used = other.used;
    allocated = other.used;
    
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
    
    return *this;
  }
  
  uninitialized_vector& operator=(uninitialized_vector&& other)
  {
    used = other.used;
    allocated = other.allocated;
    data = other.data;
    
    other.used = 0;
    other.allocated = 0;
    
    return *this;
  }
  
  uninitialized_vector (const uninitialized_vector& other) : used(other.used), allocated(other.used)
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
  
  inline size_t size() const { return used; }
  
private:
  size_t used;
  size_t allocated;
  item_t* data;
};

class own_dynamic_bitset
{
public:
  typedef unsigned int int_t;
  typedef uninitialized_vector<int_t> vec_t;
  
  static constexpr size_t sizeof_int = sizeof(int_t);
  
  own_dynamic_bitset() = default;
  ~own_dynamic_bitset() = default;
  
  inline own_dynamic_bitset(size_t size) : used(size), data(used / sizeof_int + 1) { }
  
  own_dynamic_bitset(const own_dynamic_bitset& other, size_t begin, size_t end) : used(end - begin), data(used / sizeof_int + 1)
  {
    for (int pos = 0; pos < used; pos++) {
      set(pos, other.test(pos + begin));
    }
  }
  
  inline size_t size() const { return used; }
  
  template<typename function_t>
  void generate(size_t size, function_t generator)
  {
    used = size;
    data.redim_discard(size / sizeof_int + 1);
    
    // FIXME to improve
    for (int pos = 0; pos < size; pos++) {
      set(pos, generator());
    }
  }
  
  inline void set_all(bool bit)
  {
    int_t val = 0;
    if (bit) val = ~val;
    
    for (int pos = 0; pos < data.size(); pos++)
    {
      data[pos] = val;
    }
  }
  
  inline bool test(size_t pos) const
  {
    size_t index = pos / sizeof_int;
    int subindex = pos % sizeof_int;
    int_t mask = 1 << subindex;
    
    return (data[index] & mask != 0);
  }
  inline void set(size_t pos)
  {
    size_t index = pos / sizeof_int;
    int subindex = pos % sizeof_int;
    int_t mask = 1 << subindex;
    
    data[index] |= mask;
  }
  inline void set(size_t pos, bool value)
  {
    size_t index = pos / sizeof_int;
    int subindex = pos % sizeof_int;
    int_t mask = 1 << subindex;
    
    if (value) data[index] |= mask;
    else data[index] &= ~mask;
  }
  inline void reset(size_t pos)
  {
    size_t index = pos / sizeof_int;
    int subindex = pos % sizeof_int;
    int_t mask = 1 << subindex;
    
    data[index] &= ~mask;
  }
  inline void flip(size_t pos)
  {
    size_t index = pos / sizeof_int;
    int subindex = pos % sizeof_int;
    int_t mask = 1 << subindex;
    
    data[index] ^= mask;
  }
  
  inline bool operator[](size_t pos) const
  {
    return test(pos);
  }
  
  /*
  void range_erase(size_t begin, size_t end)
  {
    //FIXME to improve
    
    own_dynamic_bitset result(used - (end - begin));
    
    int pos_dest = 0;
    for (int pos = 0; pos < used; pos++)
    {
      if (pos == begin) { pos = end - 1; continue; }
      result.set(pos_dest++, test(pos));
    }
    
    std::swap(*this, result);
  }
  
  void range_insert(size_t pos_insert, const own_dynamic_bitset& other)
  {
    //FIXME to improve
    
    own_dynamic_bitset result(used + other.used);
    
    int pos_dest = 0;
    for (int pos = 0; ; pos++)
    {
      if (pos == pos_insert) {
        for (int pos_other = 0; pos_other < other.used; pos_other++) {
          result.set(pos_dest++, other.test(pos_other));
        }
      }
      if (pos == other.used) break;
      result.set(pos_dest++, test(pos));
    }
    
    std::swap(*this, result);
  }
  */
  
  void import_string(const char* bits, size_t size)
  {
    used = size;
    data.redim_discard(size / sizeof_int + 1);
    
    //FIXME to improve
    for (size_t pos = 0; pos < size; pos++)
    {
      set(pos, (*bits++ != '0'));
    } 
  }
  
  std::string export_string() const
  {
    std::string result;
    char* buffer = new char[used + 1];
    
    //FIXME to improve
    for (size_t pos = 0; pos < used; pos++)
    {
      buffer[pos] = test(pos) ? '1' : '0';
    }
    
    buffer[used] = '\0';
    result = buffer;
    delete[] buffer;
    return result;
  }
  
private:
  size_t used;
  vec_t data;
};

