#include "bitset_safe.h"
#include "bitset_fast.h"


class own_dynamic_bitset
{
public:
  own_dynamic_bitset() = default;
  ~own_dynamic_bitset() = default;

  inline own_dynamic_bitset(size_t size)
  : bitset_safe(size)
  , bitset_fast(size)
  {
    //check("own_dynamic_bitset(size_t)");
  }

  own_dynamic_bitset(const own_dynamic_bitset& other, size_t begin, size_t end)
  : bitset_safe(other.bitset_safe, begin, end)
  , bitset_fast(other.bitset_fast, begin, end)
  {
    check("own_dynamic_bitset(own_dynamic_bitset&, begin, end)");
  }

  inline size_t size() const
  {
    size_t res = bitset_safe.size();
    if (res != bitset_fast.size()) throw "size()";
    return res;
  }

  template<typename function_t>
  own_dynamic_bitset(size_t size, function_t generator)
  : bitset_safe(size, generator)
  , bitset_fast(size)
  {
    for (size_t pos = 0; pos < size; pos++)
    {
      bitset_fast.set(pos, bitset_safe.test(pos));
    }
    
    check("generate(size_t, function_t)");
  }

  inline void set_all(bool bit)
  {
    bitset_safe.set_all(bit);
    bitset_fast.set_all(bit);
    
    check("set_all(bool)");
  }

  inline bool test(size_t pos) const
  {
    bool res = bitset_safe.test(pos);
    if (res != bitset_fast.test(pos)) throw "test(size_t)";
    return res;
  }
  inline void set(size_t pos)
  {
    bitset_safe.set(pos);
    bitset_fast.set(pos);
    
    check("set(size_t)");
  }
  inline void set(size_t pos, bool value)
  {
    bitset_safe.set(pos, value);
    bitset_fast.set(pos, value);
    
    check("set(size_t, bool)");
  }
  inline void reset(size_t pos)
  {
    bitset_safe.reset(pos);
    bitset_fast.reset(pos);
    
    check("reset(size_t)");
  }
  inline void flip(size_t pos)
  {
    bitset_safe.flip(pos);
    bitset_fast.flip(pos);
    
    check("flip(size_t)");
  }

  inline bool operator[](size_t pos) const
  {
    bool res = bitset_safe[pos];
    if (res != bitset_fast[pos]) throw "operator[](size_t)";
    return res;
  }

  void import_string(const char* bits, size_t size)
  {
    bitset_safe.import_string(bits, size);
    bitset_fast.import_string(bits, size);
    
    check("import_string(char*, size_t)");
  }

  std::string export_string() const
  {
    std::string res = bitset_safe.export_string();
    if (res != bitset_fast.export_string()) throw "export_string()";
    return res;
  }

private:
  own_dynamic_bitset_safe bitset_safe;
  own_dynamic_bitset_fast bitset_fast;

  void check(const char* context) const
  {
    std::string val_safe = bitset_safe.export_string();
    std::string val_fast = bitset_fast.export_string();
    if (val_safe != val_fast)
    {
      std::cout << "ASSERTION ERROR AFTER own_dynamic_bitset::" << context << std::endl;
      std::cout << "bitset_safe : " << val_safe << "\n";
      std::cout << "bitset_fast : " << val_fast << "\n";
      std::cout << std::endl;
      throw context;
    }
  }
};

