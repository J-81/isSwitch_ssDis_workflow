#ifndef PRED_DEL__HPP
#define PRED_DEL__HPP

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

#ifndef __INT__ONE
#define __INT__ONE
typedef signed char int1;
#endif

extern double ac_elong[21];

extern double unt;
extern double a_init;
extern double fs_n, fs_c;

#define K_MAX 5
void SetPoten(int pat, int k); // 0 - Max Sw-0.02/t^2
                               // 1 - Max Sw50+auc (change T only)
                               // 2 - (Average by position) Max Sw-0.02/t^2
                               // 3 - (Average by position) Max Sw50+auc (change T only)
                               // 4 - Null


#define PT_SIZE 100
struct PATTERN
  {
  int n;
  char pt[PT_SIZE+1];
  char& operator[](int i) {return(pt[i]);}
  };
extern int  PATTERN_T;  // 0 -- Do not use patterns
                        // 1 -- Use only HHHH pattern
                        // 2 -- Use all patterns
extern int gpt;

#define GPT  171
#define CH_PATTERN 4

#define E_PATTERN 50.

extern PATTERN pt_all[];

double Unt_Pot(void);

struct ACIDP
  {
  double e3d_n, e3d_c, p3d,
         elp_n, elp_c, plp, ept;
  double n3d;
  int  ind_ac;
  char ac;
  void reinit(void);
  void CalcP(void);
  };









inline void ACIDP::reinit(void)
  {
  e3d_n=0; e3d_c=0; p3d=0;
  elp_n=0; elp_c=0; plp=0;
  ept  =0;
  n3d=0;
  ind_ac=20;
  ac='x';
  }

class SEQ_AC
  {
  ACIDP *ac;
  int num, size;
  public:
  char info[101];
  int  calc_pattern;
  friend int Length(SEQ_AC &sa) {return(sa.num);}
  SEQ_AC(void) {ac=NULL; num=size=0; *info=0; calc_pattern=0;}
  ~SEQ_AC(void) {if(ac!=NULL) delete [] ac;}
  ACIDP& operator[](int i) {return(ac[i]);}
  void Init(int n);
  void ReInit(void);
  void LoadFasta(char *seq_fasta);
  void LoadFastaFile(const char *name);
  };

int Num_Ac(unsigned char a);
void Predict(SEQ_AC &sa);
void Predict_sh(SEQ_AC &sa);    // Предсказание без потенциалов (с одними шаблонами)
void Pred_Pattern(SEQ_AC &sa);
void PredictNC(SEQ_AC &sa);
void PredictCN(SEQ_AC &sa);
double E_Add(double a, double b);

#endif
