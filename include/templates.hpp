
#include <vector>
#include <iterator>
#include <cmath>
#include <map>

template <typename T>
class Div{
public:
  Div(double d): f(d){};
  T operator()(const T& t) const 
  {
    return t / f;
  };
private:
  double f;
};

template <typename T>
class Mul{
public:
  Mul(double d): f(d){};
  T operator()(const T& t) const
  {
    return t * f;
  };
private:
  double f;
};

template <typename T>
std::vector<T> Binning(std::vector<T>& B,int n)
{
  std::vector<T> t(B.size() / n,T());
  for(int i = 0;i < int(B.size() / n) * n;i++)
    {
      if(!(i % n))
	t[i / n] = B[i];
      else
	t[i / n] = t[i / n] + B[i];
    }
  transform(t.begin(),t.end(),t.begin(),Div<T>(n));
  return t;
}

template <typename T>
std::vector<T> Binning(std::vector<T>& B,T m,int n)
{
  std::vector<T> t(B.size() / n,T());
  for(int i = 0;i < int(B.size() / n) * n;i++)
    {
      if(!(i % n))
	t[i / n] = B[i];
      else
	t[i / n] = t[i / n] + B[i];
    }
  for(auto &i : t)
    i = i * (B.size()-1.) / (B.size() - n) - m * (n - 1.) * (B.size()) / (B.size() - n);
  return t;
}


template <typename T>
class Jackknife{
public:
  Jackknife(int n,T m): N((double)n),M(m) {};
  T operator()(const T& t) const
  { 
    if(N == 1)
	return t;
    return (M * N - t)* (1./(N - 1.));
  };
private:
  double N;
  T M;
};


template <typename T>
std::vector<T> operator*(const std::vector<T>& obj,double m){
  std::vector<T> t(obj.size(),T());
  transform(obj.begin(),obj.end(),t.begin(),Mul<T>(m));
  return t;
}


template <typename T>
std::vector<T> operator/(const std::vector<T>& obj,double m){
  std::vector<T> t(obj.size(),T());
  transform(obj.begin(),obj.end(),t.begin(),Div<T>(m));
  return t;
}

template <typename T>
class Div_Prop{
public:
  Div_Prop(int index,std::vector<std::vector<T> > prop):i(index),n(-1),p(prop) {};
  std::vector<T> operator()(std::vector<T> t)
  {
    n++;
    return t / p.at(n).at(i);
  }
private:
  int i;
  int n;
  std::vector<std::vector<T> > p;
};

template <typename T>
class SQRT{
public:
  SQRT(){};
  std::vector<T> operator()(const std::vector<T>& t)
  {
    std::vector<T> tmp;
    tmp.reserve(t.size());
    std::vector<T> s = t;
    for(typename std::vector<T>::iterator itr = s.begin();itr != s.end();itr++){
      tmp.push_back(sqrt(*itr));
    }
    return tmp;
  }
};

template <typename T>
std::vector<T> sqrt(const std::vector<T>& t)
{
  std::vector<T> temp;
  temp.reserve(t.size());

  for(T & i : t)
    temp.push_back(sqrt(i));

  return temp;
}

template <typename T>
std::vector<T> real(const std::vector<T>& t)
{
  std::vector<T> temp;
  temp.reserve(t.size());

  for(const T & i : t)
    temp.push_back(i.real());

  return temp;
}

template <typename T,typename U>
std::map<T,U> sqrt(const std::map<T,U>& t)
{
  std::map<T,U> temp;

  for(T & i : t)
    temp[i.first] = sqrt(i.second);

  return temp;
}


