#include "help.h"
#include "headers.h"
#include <chrono>

using namespace std :: chrono;
using boost::math::binomial_coefficient;
using std :: cout;
using std :: ofstream;

//--------------------------------------------------------------------------------------------------------------
//input data file paths
const string RELS_FILE = "rels_file_cpp.txt";
const string ENTS_FILE = "ents_file_cpp.txt";
const string EVIDENCE_FILE = "evidence_file_cpp.txt";
const string PROBS_FILE = "probs.txt";
const string DIM_FILE = "dims.txt";
//--------------------------------------------------------------------------------------------------------------
constexpr Index N_ITR = 50000;
constexpr Index BURN_IN = 10000;
//--------------------------------------------------------------------------------------------------------------

const Relation rels (RELS_FILE, DIM_FILE);
const Prob     prb  (PROBS_FILE);

const double Z_p[12] =
  {
    1 - 2 * prb.beta , prb.alpha         , prb.beta         , 1.0 / 3,
    prb.beta         , 1 - 2 * prb.alpha , prb.beta         , 1.0 / 3,
    prb.beta         , prb.alpha         , 1 - 2 * prb.beta , 1.0 / 3
  };

// this holds the computation of comp_prb_H and is used in comp_ch_p
// not defining it globally here causes memory fragmentation
// and slows the code dramatically
double H_p[4];

//--------------------------------------------------------------------------------------------------------------

// to be initialized by init()
Entry    ents;
vector<Index> X;


//--------------------------------------------------------------------------------------------------------------

//initialize hash tables
Hashtab relsBySrcUid = rels.group_by(R_SRCUID);
Hashtab relsByTrgUid = rels.group_by(R_TRGUID);
Hashtab relsByAppId  = rels.group_by(R_APPID);
Hashtab relsByMeshId = rels.group_by(R_MESHID);

//--------------------------------------------------------------------------------------------------------------

// this is to disambiguate the call to the overloaded log function
double (*dlog)(double) = std :: log;
int sgn (int x) { return ((x > 0) - (x < 0));}
//--------------------------------------------------------------------------------------------------------------

template<class V> void print(const V& v, std :: ostream& o=std :: cout)
{
  for (const auto& x : v)
    o << x <<"\n";
  // o << "\n";
}

std :: ostream& operator << (std :: ostream& o, const Matrix2d& m)
{
  for (size_t i=0; i<m.size1(); ++i)
    {
      for (size_t j=0; j<m.size2(); ++j)
	o << m(i,j) << " ";

      o << "\n";
    }
  return o;
}
  
//--------------------------------------------------------------------------------------------------------------

void init()
{
  // Initialize ents
  
  Matrix2s entsmat;
  init_matrix(ENTS_FILE, entsmat);

  // ents.uid = vector<Index>(ents.size1());
  // Col1s id(ents, E_UID);
  // transform(id, ents.uid, istoi);

  Col1s ents_uid(entsmat, E_UID);
  for (int j=0; j<ents_uid.size(); ++j)
    ents.uid.push_back(std::stoi(ents_uid(j)));
  
  Col1s ents_names(entsmat, E_NAME);
  auto non0_mesh = [](const string& s){ return s == "nonzeroMeSH";};
  auto n0m = find_idx_if(ents_names, non0_mesh);

  
  Col1s ents_types(entsmat, E_TYPE);
  auto mesh = [](const string& s){ return s == "MeSH";};
  auto m = find_idx_if(ents_types, mesh);
  ents.m_idx = find_difference(m, n0m);

  auto pc = [](const string& s){ return s == "Protein" || s == "Compound";};
  ents.pc_idx = find_idx_if(ents_types, pc);

  auto apl = [](const string& s){ return s == "Applicability";};
  ents.a_idx = find_idx_if(ents_types, apl);

  // read the evidence file to initialize the vector X
  Matrix2s evid;
  init_matrix(EVIDENCE_FILE, evid);
  Col1s evid_uid(evid, 0);
  Col1s evid_type(evid, 1);
  
  X = vector<Index>(ents.uid.size());
  
  auto evid_idx = match(ents_uid, evid_uid);

  for (Index i=0; i<evid_idx.size(); ++i)
    X[evid_idx[i]] = stoi(evid_type[i]);

  for (auto x : n0m)
    X[x] = 1;


}

//--------------------------------------------------------------------------------------------------------------
vector<Index> applicable_edges (const vector<Index>& net)
{
  vector<Index> res;
  
  for (auto i : net)
    if ( X[rels.rels[i+R_APPIND]] == 1 )
      res.push_back(i);

  return res;
}
//--------------------------------------------------------------------------------------------------------------
template<int N, class F, class G> array<double, N> compute (const vector<Index>& pa_idx,
							    Index ch, G g, F f)
{
  // requires f to be a functur with double operator ()(const vector<Index>&) const
  // requires g to be a functur with int operator () (int) const

  vector<Index> Val;
  for (Index i : pa_idx) Val.push_back(X[i]);
      
  vector<Index> curr_parent_idx = find_idx_if(pa_idx, equals(ch));
  array<double,N> ch_p;
  
  for (int k=0; k<N; ++k)
    {
      for (auto x : curr_parent_idx)
	Val[x] = g(k);
	  
      ch_p[k] = f(Val);
    }
  
  return apply(ch_p, dlog);
}
    
//--------------------------------------------------------------------------------------------------------------

class pc_compute
{
  const Index ch_val;
  const set<Index>& pa_idx_lab;

public:
  pc_compute(Index val, const set<Index>& s):
    pa_idx_lab{s}, ch_val{val} {}

  double operator() (const vector<Index>& val) const
  {
    vector<Index> pa_lab(pa_idx_lab.size());
    transform(pa_idx_lab, pa_lab, sgn);

    auto nmp = comp_nm_np(pa_lab, val);
    return comp_ch_p(nmp.first, nmp.second, prb.pa, ch_val);
  }
};

//--------------------------------------------------------------------------------------------------------------

void update_pc_nodes(Index t)
{
  static const Index pc_sz = ents.pc_idx.size();

  static const array<double,3> i_pc_prior_probs = {log( (1-prb.hz) / 2.0 ), log( prb.hz ), log( (1-prb.hz) / 2.0 )};
  array<double,3> i_pc_log_ch_probs = {0.0, 0.0, 0.0};

  Index i_pc = ents.pc_idx[ t % pc_sz];
  set<Index> i_pc_ch = rels.markov_blanket(relsBySrcUid.find(ents.uid[i_pc]), R_TRGUID);
    
  for (Index ch : i_pc_ch)
    {
	
      const vector<Index>& ch_net = relsByTrgUid.find(ch);
      Index ch_val = X[ rels.current_child_trg(ch_net[0]) ]; 

      vector<Index> ch_net_app = applicable_edges(ch_net);

      if ( ch_net_app.size() == 0 )
	{
	  i_pc_log_ch_probs += log(comp_ch_p(0, 0, prb.pa, ch_val));
	  continue;
	}

      if (!is_in(ents.uid[i_pc], rels.markov_blanket(ch_net_app, R_SRCUID)))
	{
	  auto nmp = comp_nm_np(ch_net_app);
	  i_pc_log_ch_probs += log(comp_ch_p(nmp.first, nmp.second, prb.pa, ch_val));
	}

      else
	{
	  set<Index> pa_idx_lab = rels.markov_blanket(ch_net_app, R_TYPE_SRCIND);
	  vector<Index> pa_idx(pa_idx_lab.size());
	  transform(pa_idx_lab, pa_idx, labs);
	  // 3 new probabilites to be computed
	  i_pc_log_ch_probs += compute<3>(pa_idx, i_pc, [](int k){return k-1;},
					  pc_compute(ch_val, pa_idx_lab));
	}

    } // end for ch
	
  array<double,3> i_pc_logX = i_pc_prior_probs + i_pc_log_ch_probs;
  
  std::initializer_list<double> i_pc_probs =
    {
      1 / (1 + exp(i_pc_logX[1] - i_pc_logX[0]) + exp(i_pc_logX[2] - i_pc_logX[0])),
      1 / (1 + exp(i_pc_logX[0] - i_pc_logX[1]) + exp(i_pc_logX[2] - i_pc_logX[1])),
      1 / (1 + exp(i_pc_logX[0] - i_pc_logX[2]) + exp(i_pc_logX[1] - i_pc_logX[2]))
    };

  //123 is the seed
  static sample i_pc_samp{123, -1};
  X[i_pc] = i_pc_samp(i_pc_probs);
}
//-----------------------------------------------------------------------------------------------------------------

void update_a_nodes(Index t){

  static const Index a_sz = ents.a_idx.size();

  Index i_a = ents.a_idx[t % a_sz];

  const vector<Index>& curr_app_net = relsByAppId.find(ents.uid[i_a]);

  vector<Index> i_a_pa_idx =  rels.get_idx(curr_app_net, R_MESHIND);
  
  //Val = X[i_a_pa_idx]
  vector<Index> Val(i_a_pa_idx.size());
  for (Index i=0; i<i_a_pa_idx.size(); ++i)
    {
      Val[i] = X[ i_a_pa_idx[i] ];
    }

  double z = 1 - pow( (1 - prb.w), count_if(Val, not_equals(0)) );
  array<double,2> i_a_pa_probs = {1-z, z};
  
  //child network. There's only one child
  set<Index> i_a_ch = rels.markov_blanket(curr_app_net, R_TRGUID);

  const vector<Index>& i_a_ch_net = relsByTrgUid.find(*i_a_ch.begin());
  
  //value of child
  Index i_a_ch_val = X[ rels.current_child_trg( i_a_ch_net[0] ) ];
  
  // temporarily set the current applicability to 1 so it would be included in the network
  X[i_a] = 1;
  
  // find the applicable edges
  vector<Index> i_a_ch_net_with_a = applicable_edges(i_a_ch_net);
  
  // temporarily set the current applicability to 0 so it would be removed from the network
  X[i_a] = 0;
  
  //find the applicable edges
  vector<Index> i_a_ch_net_without_a = applicable_edges(i_a_ch_net);
  
  // in each case compute Pr(Z | pa(Z))
  // probablility of child for different values of A
  array<double,2> i_a_ch_p;

  // set i_a_ch_p[0]
  // when removing A, are there any applicable edges?
  if ( i_a_ch_net_without_a.size() )
    {
      auto nmp = comp_nm_np(i_a_ch_net_without_a);
      i_a_ch_p[0] = comp_ch_p(nmp.first, nmp.second, prb.pc, i_a_ch_val);
    }

  else
    {
      i_a_ch_p[0] = comp_ch_p(0, 0, prb.pc, i_a_ch_val);
    }

  // set i_a_ch_p[1]
  //unique applicable network with A
  auto nmp = comp_nm_np(i_a_ch_net_with_a);
  i_a_ch_p[1] = comp_ch_p(nmp.first, nmp.second, prb.pc, i_a_ch_val);


  // computint p(A=.) = pa.probs * ch.probs / sum(pa.probs * ch.probs)
  double weight = i_a_pa_probs[0]*i_a_ch_p[0] + i_a_pa_probs[1]*i_a_ch_p[1];

  std :: initializer_list<double> i_a_probs =
    {
      i_a_pa_probs[0]*i_a_ch_p[0] / weight,
      i_a_pa_probs[1]*i_a_ch_p[1]/ weight
    };

  static sample i_a_samp{123};
  X[i_a] = i_a_samp(i_a_probs);

}

//-----------------------------------------------------------------------------------------------------------------

void update_ctx_nodes(Index t)
{

  static const Index m_sz = ents.m_idx.size();  
  Index i_m = ents.m_idx[t % m_sz];

  static const array<double,2> i_m_prior_probs = {log(prb.cp), log(1 - prb.cp)};
  
  array<double,2> i_m_log_ch_probs = {0.0, 0.0};

  set<Index> i_m_ch = rels.markov_blanket(relsByMeshId.find(ents.uid[i_m]), R_APPID);


  for (Index ch : i_m_ch)
    {
      const vector<Index>& ch_net = relsByAppId.find(ch);

      Index ch_val = X[ rels.current_child_app(ch_net[0]) ];

      vector<Index> pa_idx = rels.get_idx(ch_net, R_MESHIND);
      // update 2 new probabilities
      i_m_log_ch_probs += compute<2>(pa_idx, i_m, [](int i){return i;},
				     [ch_val](const vector<Index>& Val)
				     {
				       double z = 1 - pow((1.0 - prb.w), count_if(Val, not_equals(0)));
				       return ch_val == 1 ? z : 1 - z;
				     });				     
    }
    
  array<double,2> i_m_logX = i_m_prior_probs + i_m_log_ch_probs;
    
  std::initializer_list<double> i_m_probs =
    {
      1 / (1 + exp(i_m_logX[1] - i_m_logX[0])),
      1 / (1 + exp(i_m_logX[0] - i_m_logX[1]))
    };
  
  static sample i_m_samp{123};
  X[i_m] = i_m_samp(i_m_probs);
  
}// update_ctx_nodes

//----------------------------------------------------------------------------------------------------------------

class ma_marginals
{
  
protected:
  Matrix2d m;
  const vector<Index>& v;
  
public:
  ma_marginals (Index d1, const vector<Index>& idx) :
    v{idx}
  {
    m = Zmatrix2d(d1, idx.size());
  }
  
  void update()
  {
    for (Index i=0; i<m.size2(); ++i)
      X[v[i]] == 0 ? ++m(0,i) : ++m(1,i);
  }

  void write (const string& file)
  {
    std::ofstream fout(file);
    Matrix2d tr = trans(m);
    tr *= 1.0 / N_ITR;
    fout << tr;
    fout.close();
  }
};

class pc_marginals : public ma_marginals
{
public:

  pc_marginals (Index d1, const vector<Index>& idx) :
    ma_marginals(d1,idx) {}

  void update()
  {
    for (Index i=0; i<m.size2(); ++i)

      if (X[v[i]] == 0)
	++m(1,i);
      else if (X[v[i]] == 1)
	++m(2,i);
      else
	++m(0,i);
  }   
};

//-----------------------------------------------------------------------------------------------------------------

int main()
{


  try
    {
      auto t1 = steady_clock :: now();
      init();
      
      const Index e1 = BURN_IN + 1;
      const Index e2 = BURN_IN + N_ITR + 1;

      std::cout << "Before Burn in ...\n";
      for(Index t=1; t != e1; ++t)
	{

	  if (t%50000 == 0)
	    {
	      std::cout <<"Iteration: " << t << " / " << N_ITR <<"\n";
	    }

	  update_pc_nodes (t);
	  update_a_nodes  (t);
	  update_ctx_nodes(t);
	}

      std::cout <<"After Burn in ...\n";
      // Marginal probabilities will be updated.
      
      ma_marginals marg_m(2,ents.m_idx);
      ma_marginals marg_a(2,ents.a_idx);
      pc_marginals marg_pc(3,ents.pc_idx);
      
      for(Index t=e1; t != e2; ++t)
	{
	  if (t%100000 == 0)
	    {
	      std::cout <<"Iteration: " << t << " / " << N_ITR <<"\n";
	    }
	  
	  update_pc_nodes (t);
	  marg_pc.update();
	  
	  update_a_nodes  (t);
	  marg_a.update();
	  
	  update_ctx_nodes(t);
	  marg_m.update();
	  
	}

      auto t2 = steady_clock::now() - t1;
      std::cout << "\nSimulation Time:"
		  <<  duration_cast<seconds>(t2).count()
		  << "(s)\n";

      marg_pc.write("marg_prob_pc.txt");
      marg_m.write("marg_prob_m.txt");
      marg_a.write("marg_prob_a.txt");
  
    }
  
  catch(Hash_error h)
    {
      std:: cout << h.name;
    }
  
  return 0;
}

std::pair<double, double> comb(int np, int nm, int k, int l, double pc){

    std::pair<double,double> res;

    int n = np + nm;

    double x = binomial_coefficient<double>(np, k) *
      binomial_coefficient<double>(nm, l) *
      pow(1 - pc, k+l) *
      pow(pc, n-k-l);
    double y = static_cast<double>(abs(k-l))/((k+l)*(k+l));
    res.first = x * k * y;
    res.second = x * l * y;
    return res;
  }

/*------------------------------------------------------------------------------------------
 * This function computes P(H | X)
 */

void comp_prob_H(int nm, int np, double c){

    int n = nm + np;
    
    if ( n <= 0 )
      {
	H_p[0] = prb.pm;
	H_p[1] = prb.pz;
	H_p[2] = prb.pp;
	H_p[3] = 0;
	return;
      }

    double pcn = pow(c, n);
    double phm = pcn * prb.pm;
    double ph0 = pcn * prb.pz;
    double php = pcn * prb.pp;

    if ( nm >= 1 )
      phm += pow(c, np) - pcn;
    
    if ( np >= 1 )
      php += pow(c, nm) - pcn;
    
    double B = 0.0;
    double C = 0.0;
    
    if ( np >= 1 && nm >= 1 )
      {
	for (int k=1; k<=np; ++k)
	  {
	    for(int l=1; l<=nm; ++l)
	      {
		auto x = comb(np, nm, k, l, c);
		B += x.first;
		C += x.second;
	      }
	  }
      }

    
    phm += C;
    php += B;
    double pha = (1 - ( pow(c, nm) + pow(c, np) - pcn )) - B - C;

    H_p[0] = phm;
    H_p[1] = ph0;
    H_p[2] = php;
    H_p[3] = pha;
}
  

/*------------------------------------------------------------------------------------------
 * This function computes the child probabilies using compProbH and Z_p
 */

double comp_ch_p(int nm, int np, double pac, int val){
 
    comp_prob_H(nm, np, pac * (1-abs(val)) + prb.pc * abs(val));
    //return  inner_prod(Z_p, H_p);
    int i = ++val << 2;
    return
      H_p[0] * Z_p[i] + H_p[1] * Z_p[i+1] + H_p[2] * Z_p[i+2] + H_p[3] * Z_p[i+3];
      
  }

/*------------------------------------------------------------------------------------------*/
  
std :: pair<Index, Index> comp_nm_np(const vector<Index>& Val, const vector<Index>& pa_lab)
{

  Index sz= Val.size();
  vector<Index> pred_vals(sz);

  for (Index i=0; i<sz; ++i)
    pred_vals[i] = Val[i] * pa_lab[i];

  std::pair<Index,Index> res(0,0);
  res.first  = count_if(pred_vals, equals(-1)); //nm
  res.second = count_if(pred_vals, equals(1)); //np

  return res;
}
/*------------------------------------------------------------------------------------------*/

std::pair<Index, Index> comp_nm_np(const vector<Index>& n)
{

  const set<Index> pa_ind_lab = rels.markov_blanket(n, R_TYPE_SRCIND);
  auto sz = pa_ind_lab.size();

  vector<Index> pa_ind(sz);
  vector<Index> pa_lab(sz);
  vector<Index> Val(sz);

  transform(pa_ind_lab, pa_ind, labs);
  transform(pa_ind_lab, pa_lab, sgn);
  

  for (Index i=0; i<sz; ++i) Val[i] = X[pa_ind[i]];

  return comp_nm_np(Val, pa_lab);
}
/*------------------------------------------------------------------------------------------*/
