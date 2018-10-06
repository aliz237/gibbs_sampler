#include "headers.h"

//----------------------------------------------------------------------------------------------------------
void Hashtab :: insert (Index key, Index val)
{
  Index i = hash_map (key);

  auto pos = std :: find_if (tab[i].begin(), tab[i].end(),
			     [key](const Node& it){ return it.key == key; });
  
  if (pos != tab[i].end())
    pos->idx.push_back(val);
  
  else
    {
      Node n;
      n.idx.push_back(val);
      n.key = key;
      tab[i].push_back(n);
    }
}
//----------------------------------------------------------------------------------------------------------

vector<Index> Hashtab :: find (Index key)
{
  Index i = hash_map (key);

  auto pos = find_if (tab[i].begin(), tab[i].end(),
		     [key](const Node& n){ return n.key == key; });
  
  if (pos == tab[i].end())
    herror("can't find" + std :: to_string(key) + "in hash table.\n");
  
  return pos->idx;
}

//----------------------------------------------------------------------------------------------------------

Index Hashtab :: next_prime(Index i)
{
  if (i < 2) return 2;
  if (i % 2 == 0) ++i;

  while (1)
    {
      if (is_prime(i))
	return i;
      
      i += 2;
    }
  /*
  for (; !is_prime(i); i += 2;);
  return i;
  */
}

//----------------------------------------------------------------------------------------------------------

bool Hashtab :: is_prime (Index i)
{
  // this is called with i > 2
  Index end = ceil(sqrt(i));
  
  for (Index j=3; j<=end; j += 2)
    if ( i % j == 0)
      return false;
  
  return true;
}    

//----------------------------------------------------------------------------------------------------------
unsigned int Hashtab :: hash_map2( unsigned int a)
{
   a = (a+0x7ed55d16) + (a<<12);
   a = (a^0xc761c23c) ^ (a>>19);
   a = (a+0x165667b1) + (a<<5);
   a = (a+0xd3a2646c) ^ (a<<9);
   a = (a+0xfd7046c5) + (a<<3);
   a = (a^0xb55a4f09) ^ (a>>16);
   return a % tab.size();
}
