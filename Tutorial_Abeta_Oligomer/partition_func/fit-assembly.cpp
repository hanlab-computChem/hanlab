#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<vector>
#include<string>
#include<math.h>

using std::vector;
using std::cout;


#define  show_array(___a, ___b)  \
for(int ___i=0;___i<(___a);___i++) \
 cout<< (___b)[___i]  <<" "; \
 cout<<"\n";


double diff_ns(int n_cut, vector<double> & _ns_ref, 
			vector<double> & _ns_mean ){

 double _r=0.0;

 for (int _i=0; _i<n_cut; _i++)
  _r+=(log(_ns_ref[_i]/_ns_mean[_i])*log(_ns_ref[_i]/_ns_mean[_i]));

 return _r/(double)n_cut;
}


//long int _par_n(int n, vector<vector<int > > & _par_s,
//			vector<vector<int> > & _par_ns){

long int _par_n(int n, int n_cut,  vector<double> & _q0,
			vector<double> & _ns_mean){




if (_q0.size()!= n) {cout<<"Error! the number of inital guess of q0 is not "<<n<<".\n";exit(1);}

vector<int> _x(n,1);
_ns_mean.clear();
_ns_mean.resize(n,0.0);
//these are predetermind values for efficiency
vector<double> _nk_r(n);
vector<int> _ns(n);

for(int _i=0;_i<n; _i++) _nk_r[_i] = 1.0/(double)(_i + 1);


 
 _x[0]= n;
 int m = 1;
 int h =1 ;
 int r,t;
 long int cnt=1;
 double _Q=0.0;

// for(int _i=0;_i<m;_i++)
//   cout<<  _x[_i]<<" ";
// cout<<"\n";

bool _to_eval;

 _to_eval = true;

for (int _i=0; _i<m; _i++)
 if (_x[_i]>n_cut) {_to_eval=false; break;}
if (_to_eval)
{
  int _xi;
  int __t = -1;
  int __c = 0;
  for(int _i=0;_i<n;_i++) _ns[_i]=0;

  double _w=1.0;  

for(int _i=0;_i<m;_i++)
{
  _xi = _x[_i];
  if (_xi!=__t) {__c=0;__t=_xi;}
  __c++;
  _w*=(_q0[_xi-1] * _nk_r[__c-1]);
  _ns[_xi-1]++; 
//  cout<< _xi  <<""<< __c <<" " << _q0[_xi-1] <<" "<< _nk_r[__c-1]<<"\n";
}
// cout<<"\n";
 _Q+=_w;

for(int _i=0;_i<n;_i++)
 _ns_mean[_i]+=((double)_ns[_i] * _w);


}
// vector<int> __s,__ns;
/* 
 {
  int __t = -1;
  int __c = 0;

 __t = _x[0];
 __s.push_back(__t);
 __c++;

 for(int _i=1;_i<m;_i++)

  {
    if (_x[_i]!=__t) {
      __ns.push_back(__c);
      __t = _x[_i];
      __s.push_back(__t); 
      __c = 0; }
    __c++;
    
  } 
  __ns.push_back(__c);

 _par_s.push_back(__s);
 _par_ns.push_back(__ns);

 }
 
*/

 while (_x[0]!=1){
  if (_x[h-1]==2) {
   m+=1;
   _x[h-1]=1;
   h-=1; }
  else {
   r= _x[h-1]-1;
   t = m - h +1;
   _x[h-1]=r;
   while (t>=r) {
    h+=1;
    _x[h-1]=r;
    t-=r; }
   if (t==0){
    m=h; }
   else {
    m=h+1;
    if (t>1) {
     h+=1;
     _x[h-1]=t;}
    }
   }
  cnt++;

  if (cnt%1000000==0) cout<<cnt/1000000<<"M combs searched..\n"<< std::flush;
//  for(int _i=0;_i<m;_i++) 
//   cout<<  _x[_i]<<" ";
// cout<<"\n";
/*
 {
  int __t = -1;
  int __c = 0;
 __s.clear();
 __ns.clear(); 
 __t = _x[0];
 __s.push_back(__t);
 __c++;
 
 for(int _i=1;_i<m;_i++)
  
  { 
    if (_x[_i]!=__t) {
      __ns.push_back(__c);
      __t = _x[_i];
      __s.push_back(__t);
      __c = 0; }
    __c++;
   
  } 
  __ns.push_back(__c);

 _par_s.push_back(__s);
 _par_ns.push_back(__ns);

 }
*/
 _to_eval = true;

for (int _i=0; _i<m; _i++)
 if (_x[_i]>n_cut) {_to_eval=false; break;}
if (_to_eval)
{
  int _xi;
  int __t = -1;
  int __c = 0;
  for(int _i=0;_i<n;_i++) _ns[_i]=0;

  double _w=1.0;

for(int _i=0;_i<m;_i++)
{
  _xi = _x[_i];
  if (_xi!=__t) {__c=0; __t = _xi;}
  __c++;
  _w*=(_q0[_xi-1] * _nk_r[__c-1]);
  _ns[_xi-1]++;

//  cout<< _xi <<" "<< __c <<" " << _q0[_xi-1] <<" "<< _nk_r[__c-1]<<" "<< _w <<"\n";
}

 _Q+=_w;

// cout<< _Q <<"\n";

for(int _i=0;_i<n;_i++)
 _ns_mean[_i]+=((double)_ns[_i] * _w);


}
/*
for(int _i=0;_i<n;_i++) cout<< _ns[_i] <<" ";
cout<<"\n";
for(int _i=0;_i<n;_i++) cout<< _ns_mean[_i] <<" ";
cout<<"\n";
cout<<"\n";
*/

 }

 for (int _i=0;_i<n;_i++) _ns_mean[_i]/=_Q;

 return cnt; 

}



int main(int argc, const char ** argv){


int n =100;
int n_cut = 28;
vector<double> _q0(n,1.0); 
vector<double> _ns_mean(n,1.0);
vector<double> _ns_ref(n,1.0);
double _e = 1e-6;
int _cycle =0;

/*
_ns_ref[0]= 0.825415;
_ns_ref[1]= 1.35878;
_ns_ref[2]= 0.553996;
_ns_ref[3]= 0.449426;
_ns_ref[4]= 0.362568;
_ns_ref[5]= 0.290727;
_ns_ref[6]= 0.231576;
_ns_ref[7]= 0.183117;
_ns_ref[8]= 0.143622;
_ns_ref[9]= 0.111643;
_ns_ref[10]= 0.085871;
_ns_ref[11]= 0.0653577;
_ns_ref[12]= 0.0489319;
_ns_ref[13]= 0.0363795;
_ns_ref[14]= 0.025968;
_ns_ref[15]= 0.019136;
_ns_ref[16]= 0.0120193;
_ns_ref[17]= 0.00948892;
*/


//ab40 n_cut=23

_ns_ref[	0	]=	7.0206;
_ns_ref[	1	]=	22.1080;
_ns_ref[	2	]=	26.7352;
_ns_ref[	3	]=	13.3373;
_ns_ref[	4	]=	6.7018;
_ns_ref[	5	]=	8.3774;
_ns_ref[	6	]=	3.6809;
_ns_ref[	7	]=	2.0331;
_ns_ref[	8	]=	3.8411;
_ns_ref[	9	]=	1.4837;
_ns_ref[	10	]=	1.6959;
_ns_ref[	11	]=	1.2094;
_ns_ref[	12	]=	0.4301;
_ns_ref[	13	]=	0.3410;
_ns_ref[	14	]=	0.3816;
_ns_ref[	15	]=	0.0892;
_ns_ref[	16	]=	0.1864;
_ns_ref[	17	]=	0.0871;
_ns_ref[	18	]=	0.0148;
//_ns_ref[	19	]=	0.2246;
//_ns_ref[	20	]=	0.0062;
//_ns_ref[	21	]=	0.0060;
//_ns_ref[	22	]=	0.0042;


//ab42 n_cut=28
/*
_ns_ref[	0	]=	6.4021;
_ns_ref[	1	]=	31.2574;
_ns_ref[	2	]=	20.4904;
_ns_ref[	3	]=	10.4788;
_ns_ref[	4	]=	7.1987;
_ns_ref[	5	]=	2.6164;
_ns_ref[	6	]=	0.8796;
_ns_ref[	7	]=	1.6424;
_ns_ref[	8	]=	1.1524;
_ns_ref[	9	]=	2.8304;
_ns_ref[	10	]=	2.4353;
_ns_ref[	11	]=	1.2852;
_ns_ref[	12	]=	3.3062;
_ns_ref[	13	]=	1.3984;
_ns_ref[	14	]=	2.4175;
_ns_ref[	15	]=	0.7550;
_ns_ref[	16	]=	0.4482;
_ns_ref[	17	]=	0.8388;
_ns_ref[	18	]=	0.1448;
_ns_ref[	19	]=	0.4138;
//_ns_ref[	20	]=	1.3121;
//_ns_ref[	21	]=	0.1391;
//_ns_ref[	22	]=	0.0579;
//_ns_ref[	23	]=	0.0227;
//_ns_ref[	24	]=	0.0497;
//_ns_ref[	25	]=	0.0107;
//_ns_ref[	26	]=	0.0024;
//_ns_ref[	27	]=	0.0099;
*/

int _larg_clu = sizeof(_ns_ref); 
for (int _i=1;_i<_larg_clu+1;_i++) 
 _ns_ref[_i-1]/=_i;

double _d = diff_ns(n_cut, _ns_ref, _ns_mean);
double _dq = 1.0;
while (_dq>_e){

 cout<<"cycle "<<_cycle<<"\n";

 

 cout<<"q0 : ";
 show_array(n_cut, _q0) 
 _par_n(n,n_cut, _q0, _ns_mean);

 _cycle++;
 _d = diff_ns(n_cut, _ns_ref, _ns_mean);

 cout<<"_ns_ref : ";
 show_array(n_cut, _ns_ref)
 cout<<"_ns_fit : ";
 show_array(n_cut, _ns_mean)
 cout<<"d : "<< _d <<"\n";;

 cout<<"updating q0...\n";
 _dq=0.0;
 for(int _i=1;_i<n_cut;_i++) {
  double _q_p = _q0[_i]*pow(_ns_ref[_i]/_ns_mean[_i]*pow(_ns_mean[0]/_ns_ref[0], _i+1 ), 0.75  );
  _dq+= (log(_q_p/_q0[_i])* log(_q_p/_q0[_i]) );
  _q0[_i] = _q_p;
 }
 cout<<"dq : "<< _dq  <<"\n";
 cout<<"------\n";
 //break;
}

 cout<<"Converged q0: ";
 show_array(n_cut, _q0)

}
