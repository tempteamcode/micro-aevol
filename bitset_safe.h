#include <vector>
#include <utility>

class own_dynamic_bitset_safe
{
public:
  own_dynamic_bitset_safe() = default;
  ~own_dynamic_bitset_safe() = default;

  inline own_dynamic_bitset_safe(size_t size) : used(size), data(used+1,'0') { data[used] = '\0'; }

  own_dynamic_bitset_safe(const own_dynamic_bitset_safe& other, size_t begin, size_t end) : used(end - begin), data(other.data.begin() + begin, other.data.begin() + end + 1)
  {
    data[used] = '\0';
  }

  inline size_t size() const { return used; }

  template<typename function_t>
  inline own_dynamic_bitset_safe(size_t size, function_t generator) : used(size), data(used+1,'0')
  {
    for (int pos = 0; pos < size; pos++) {
      set(pos, generator());
    }
    
    data[used] = '\0';
  }

  inline void set_all(bool bit)
  {
    char val = bit ? '1' : '0';
    
    for (int pos = 0; pos < used; pos++)
    {
      data[pos] = val;
    }
  }

  inline bool test(size_t pos) const
  {
    return data[pos] != '0';
  }
  inline void set(size_t pos)
  {
    data[pos] = '1';
  }
  inline void set(size_t pos, bool value)
  {
    data[pos] = (value ? '1' : '0');
  }
  inline void reset(size_t pos)
  {
    data[pos] = '0';
  }
  inline void flip(size_t pos)
  {
    data[pos] = (data[pos] == '0' ? '1' : '0');
  }

  inline bool operator[](size_t pos) const
  {
    return test(pos);
  }

  void import_string(const char* bits, size_t size)
  {
    used = size;
    data = std::vector<char>(bits, bits+size);
  }

  std::string export_string() const
  {
    return std::string(data.begin(), data.end()-1);
  }

private:
  size_t used;
  std::vector<char> data;
};

