
#ifndef SORT_HPP
#define SORT_HPP

#include <vector>

template<typename T>
std::vector<T> forward(std::vector<T> t,int s)
{
  std::vector<T> temp;
  
  for(int i = 0;i < t.size();i++)
    temp.push_back(t.at((i + s) % t.size()));
  return temp;
}

template<typename T>
std::vector<T> backward(std::vector<T> t,int s)
{
  std::vector<T> temp;
  
  for(int i = 0;i < t.size();i++)
    temp.push_back(t.at((2 * t.size() - i + s) % t.size()));
  return temp;
}


#endif
