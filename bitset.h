#ifndef OWN_BITSET_H
#define OWN_BITSET_H

#include <vector>
#include <utility>
#include <cstring> // std::memcpy

#include "bits.h"

class own_bitset
{
public:
  typedef int size_t;
  typedef uint32_t int_t;
  static constexpr size_t sizeof_int = sizeof(int_t) * 8;
  
  ~own_bitset() { if (data != nullptr) delete[] data; }

  inline own_bitset(const own_bitset& other) : used(other.used), data(new int_t[(used+sizeof_int-1) / sizeof_int])
  {
    std::memcpy(data, other.data, ((used+sizeof_int-1) / sizeof_int) * sizeof(int_t));
  }

  own_bitset(own_bitset&& other) : used(other.used), data(other.data)
  {
    const_cast<int_t*&>(other.data) = nullptr;
  }

  inline own_bitset& operator=(const own_bitset& other)
  {
    if (&other == this) return *this;
    if (data != nullptr) delete[] data;
    const_cast<size_t&>(used) = other.used;
    const_cast<int_t*&>(data) = new int_t[(used+sizeof_int-1) / sizeof_int];
    std::memcpy(data, other.data, ((used+sizeof_int-1) / sizeof_int) * sizeof(int_t));
    return *this;
  }

  inline own_bitset& operator=(own_bitset&& other)
  {
    if (&other == this) return *this;
    if (data != nullptr) delete[] data;
    const_cast<size_t&>(used) = other.used;
    const_cast<int_t*&>(data) = other.data;
    //const_cast<size_t&>(other.used) = 0;
    const_cast<int_t*&>(other.data) = nullptr;
    return *this;
  }

  inline own_bitset(size_t size) : used(size), data(new int_t[(used+sizeof_int-1) / sizeof_int]) { }

  template<typename function_t>
  own_bitset(size_t size, function_t generator) : used(size), data(new int_t[(used+sizeof_int-1) / sizeof_int])
  {
    size_t count = used;
    int_t* ptr = data;
    int_t val = 0;
    int_t mask = 1;
    while (count --> 0) {
      if (generator()) val |= mask;
      if ((mask <<= 1) == 0) {
        *ptr++ = val;
        val = 0;
        mask = 1;
      }
    }
    if (mask != 1) *ptr++ = val;
  }

  own_bitset(const char* bits, size_t size) : used(size), data(new int_t[(used+sizeof_int-1) / sizeof_int])
  {
    size_t count = used;
    int_t* ptr = data;
    int_t val = 0;
    int_t mask = 1;
    while (count --> 0) {
      if (*bits++ == '1') val |= mask;
      if ((mask <<= 1) == 0) {
        *ptr++ = val;
        val = 0;
        mask = 1;
      }
    }
    if (mask != 1) *ptr++ = val;
  }

  inline size_t size() const { return used; }

  inline void set_all(bool bit)
  {
    size_t size = ((used+sizeof_int-1) / sizeof_int);
    int_t val = 0;
    if (bit) val = ~val;
    
    for (size_t index = 0; index < size; index++)
    {
      data[index] = val;
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

  std::string export_string() const
  {
    std::string result;
    char* buffer = new char[used + 1];
    char* bits = buffer;

    int_t* ptr = data;
    size_t count = used;
    int_t val = 0;
    int_t mask = 0;
    while (count --> 0) {
      if ((mask <<= 1) == 0) {
        val = *ptr++;
        mask = 1;
      }
      *bits++ = (val & mask) ? '1' : '0';
    }

    *bits = '\0';
    result = buffer;
    delete[] buffer;
    return result;
  }

  int_t getSequenceRev(int pos, int len) const {
      if (pos >= used) pos -= used;
      size_t index = (pos/ sizeof_int);
      int subindex = pos % sizeof_int;

      int_t res;

      if(subindex + len <= sizeof_int && (pos%used)+len <= used){
          int_t mask = lookupmasks[len-1];

          mask = mask << subindex;

          res = (mask&data[index])>>subindex;

      }else if(subindex + len > sizeof_int && pos+len > used) {
          int len1 = sizeof_int - subindex;
          int len2 = used - (index + 1) * sizeof_int;
          int len3 = len - len1 - len2;

          int_t mask1 = lookupmasks[len1-1];

          mask1 = mask1 << subindex;

          res = (mask1 & data[index]) >> (subindex);

          int_t mask2 = lookupmasks[len2-1];

          res |= (mask2 & data[index + 1]) << len1;

          int_t mask3 = lookupmasks[len3-1];

          res |= (mask3 & data[0]) << (len1 + len2);
      }else{
          int size_left;

          if(pos+len <= used)
              size_left = sizeof_int-subindex;
          else
              size_left = used - index*sizeof_int - subindex;

          int len1 = size_left;
          int len2 = len - len1;

          int_t mask1 = lookupmasks[len1-1];

          mask1 = mask1 << subindex;

          res = (mask1&data[index])>>(subindex);

          int_t mask2 = lookupmasks[len2-1];

          int newIndex = index+1;
          if(pos+len > used)
              newIndex = 0;

          res |= (mask2&data[newIndex])<<len1;
      }

      return res;
  }

  inline int_t getSequence(int pos, int len) const {
    return bits_reverse(len, getSequenceRev(pos, len));
  }

  template <typename callback_t>
  int_t forSequences(size_t start, size_t count, int len, callback_t callback) const {
    int_t result = getSequence(start, len);
    int_t mask = lookupmasks[len-1];

    size_t pos = start + len;
    if (pos >= used) pos -= used;
    callback(pos, result);

    size_t index = pos / sizeof_int;
    int subindex = pos % sizeof_int;
    int_t* ptr = &data[index];
    int_t val = *ptr >> subindex;

    if (count --> 0)
    while (count --> 0)
    {
      result = ((result << 1) & mask) | (val & 1);
      val >>= 1;

      if (++pos == used)
      {
        index = 0; subindex = 0; val = *(ptr = data);
      }
      else if (++subindex == sizeof_int)
      {
        index++; subindex = 0; val = *++ptr;
      }

      callback(pos, result);
    }
  }

private:
  const size_t used;
  int_t* const data;
};



#endif //OWN_BITSET_H
