#ifndef OWN_BITSET_H
#define OWN_BITSET_H

#include <vector>
#include <utility>
#include <cstring> // std::memcpy

#include "bits.h"

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
        std::memcpy(data, other.data, used * sizeof(item_t));
      } else {
        item_t* dest = data;
        const item_t* src = other.data;
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
    //other.data = nullptr;
    
    return *this;
  }
  
  uninitialized_vector (const uninitialized_vector& other) : used(other.used), allocated(other.used)
  {
    if (allocated > 0) {
      data = new item_t[allocated];
      
      if (std::is_trivially_copyable<item_t>::value) {
        std::memcpy(data, other.data, used * sizeof(item_t));
      } else {
        item_t* dest = data;
        const item_t* src = other.data;
        size_t count = used;
        while (count --> 0) *dest++ = *src++;
      }
    }
  }
  
  uninitialized_vector (uninitialized_vector&& other) : used(other.used), allocated(other.allocated), data(other.data)
  {
    other.used = 0;
    other.allocated = 0;
    //other.data = nullptr;
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
            std::memcpy(data, data_old, used * sizeof(item_t));
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
  typedef uint32_t int_t;
  typedef uninitialized_vector<int_t> vec_t;
  
  static constexpr size_t sizeof_int = sizeof(int_t) * 8;
  
  own_dynamic_bitset() = default;
  ~own_dynamic_bitset() = default;
  
  inline own_dynamic_bitset(size_t size) : used(size), data(used / sizeof_int + 1) { }
  
  own_dynamic_bitset(const own_dynamic_bitset& other, size_t begin, size_t end) : used(end - begin), data(used / sizeof_int + 1)
  {
    for (int pos = 0; pos < used; pos++) {
      set(pos, other.test(pos + begin));
    }
  }
  
  template<typename function_t>
  inline own_dynamic_bitset(size_t size, function_t generator) : used(size), data(used / sizeof_int + 1)
  {
    /*
    for (size_t pos = 0; pos < size; pos++) {
      set(pos, generator());
    }
    */
    
    size_t count = used;
    size_t index = 0;
    int_t val = 0;
    int_t mask = 1;
    while (count --> 0) {
      if (generator()) val |= mask;
      if ((mask <<= 1) == 0) {
        data[index++] = val;
        val = 0;
        mask = 1;
      }
    }
    if (mask != 1) data[index++] = val;
  }
  
  inline size_t size() const { return used; }
  
  inline void set_all(bool bit)
  {
    int_t val = 0;
    if (bit) val = ~val;
    
    for (size_t pos = 0; pos < data.size(); pos++)
    {
      data[pos] = val;
    }
  }
  
  inline bool test(size_t pos) const
  {
    size_t index = pos / sizeof_int;
    int subindex = pos % sizeof_int;
    int_t mask = 1 << subindex;
    
    return ((data[index] & mask) != 0);
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

  void import_string(const char* bits, size_t size)
  {
    used = size;
    data.redim_discard(size / sizeof_int + 1);
    
    /*
    for (size_t pos = 0; pos < size; pos++)
    {
      set(pos, (*bits++ != '0'));
    }
    */
    
    size_t count = used;
    size_t index = 0;
    int_t val = 0;
    int_t mask = 1;
    while (count --> 0) {
      if (*bits++ == '1') val |= mask;
      if ((mask <<= 1) == 0) {
        data[index++] = val;
        val = 0;
        mask = 1;
      }
    }
    if (mask != 1) data[index++] = val;
  }
  
  char* export_string(char* bits, size_t size) const
  {
    if (size < used + 1)
    {
      delete[] bits;
      bits = new char[used + 1];
    }
    char* buffer = bits;

    /*
    for (size_t pos = 0; pos < used; pos++)
    {
      buffer[pos] = (test(pos) ? '1' : '0');
    }
    */
    
    size_t count = used;
    size_t index = 0;
    int_t val = 0;
    int_t mask = 0;
    while (count --> 0) {
      if ((mask <<= 1) == 0) {
        val = data[index++];
        mask = 1;
      }
      *bits++ = (val & mask) ? '1' : '0';
    }
    
    *bits = '\0';
    return buffer;
  }

  inline std::string export_string() const
  {
    char* buffer = export_string(nullptr, 0);
    std::string result(buffer);
    delete[] buffer;
    return result;
  }

  int_t getSequenceRev(int pos, int len) {
      if (pos >= used) pos -= used;
      size_t index = (pos/ sizeof_int);
      int subindex = pos % sizeof_int;

      int_t res;

      if(subindex + len <= sizeof_int && (pos%used)+len <= used){
          int_t mask = MASKS[len-1];

          mask = mask << subindex;

          res = (mask&data[index])>>subindex;

      }else if(subindex + len > sizeof_int && pos+len > used) {
          int len1 = sizeof_int - subindex;
          int len2 = used - (index + 1) * sizeof_int;
          int len3 = len - len1 - len2;

          int_t mask1 = MASKS[len1-1];

          mask1 = mask1 << subindex;

          res = (mask1 & data[index]) >> (subindex);

          int_t mask2 = MASKS[len2-1];

          res |= (mask2 & data[index + 1]) << len1;

          int_t mask3 = MASKS[len3-1];

          res |= (mask3 & data[0]) << (len1 + len2);
      }else{
          int size_left;

          if(pos+len <= used)
              size_left = sizeof_int-subindex;
          else
              size_left = used - index*sizeof_int - subindex;

          int len1 = size_left;
          int len2 = len - len1;

          int_t mask1 = MASKS[len1-1];

          mask1 = mask1 << subindex;

          res = (mask1&data[index])>>(subindex);

          int_t mask2 = MASKS[len2-1];

          int newIndex = index+1;
          if(pos+len > used)
              newIndex = 0;

          res |= (mask2&data[newIndex])<<len1;
      }

      return res;
  }

  inline int_t getSequence(int pos, int len) {
    return bits_reverse(len, getSequenceRev(pos, len));
  }

private:
  size_t used;
  vec_t data;
};

#endif //OWN_BITSET_H
