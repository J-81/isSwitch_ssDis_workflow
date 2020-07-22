#include "pred_del.hpp"
                     // Max Sw50+auc (by temperature), All pattern
                     // T iter=15/  4 Sw= 0.48387703    Sw_Tpot=  0.47262703  Tpot=  0.75000000  Q2= 0.76630477 Q2_Tpot=  0.75505477
                     //     CYS,     MET,     PHE,     ILE,     LEU,     VAL,     TRP,     TYR,     ALA,     GLY,     THR,     SER,     GLN,     ASN,     GLU,     ASP,     HIS,    ARG,      LYS,    PRO,      UNK
double ac_elong[21] = {  1.0980,  0.0899,  1.0183,  0.9935,  0.6880,  0.7298,  1.5734,  1.0513,  0.0289, -0.1192, -0.0228, -0.6032, -0.4020, -0.2612, -0.4826, -0.3638, -0.0110, -0.2658, -0.3660, -0.5851,  0.0000};
double a_init= 3.4848;
double fs_n  =-4.3152;
double fs_c  =-4.2289;

const double wac_elong[21] =
  //     CYS,      MET,      PHE,      ILE,      LEU,      VAL,      TRP,      TYR,      ALA,      GLY,      THR,      SER,      GLN,      ASN,      GLU,      ASP,      HIS,     ARG,       LYS,     PRO,       UNK
  {0.0137229,0.0224214,0.0389551,0.0545193,0.0887486,0.0688628,0.0141371,0.0346998,0.0794769,0.0742804,0.0553257,0.0608747,0.0369777,0.0424005,0.0654720,0.0569914,0.0240030,0.0493192,0.0585566,0.0461315,0.0000341};

const double wa_init=0.0064994;
const double wa_fs  =0.0037950;
double unt=1;

int  PATTERN_T=2;
//int gpt=  0; // PATTERN_T == 0 -- Do not use patterns
//int gpt=  1; // PATTERN_T == 1 -- Use only HHHH pattern
int   gpt=GPT; // PATTERN_T == 2 -- Use all patterns

PATTERN pt_all[GPT] = {
   { 4, "HHHH"                 },  //   1
   { 6, "ENLYFQ"               },  //   2
   { 4, "GSHM"                 },  //   3
   { 5, "GPGSM"                },  //   4
   { 5, "DDDDK"                },  //   5
   { 8, "EGGHHHHH"             },  //   6
   { 5, "VPRGS"                },  //   7
   {11, "ASMTGGQQMGR"          },  //   8
   { 6, "LEAHHH"               },  //   9
   { 4, "SNAM"                 },  //  10
   { 7, "GGGGSGG"              },  //  11
   { 5, "GPLGS"                },  //  12
   { 8, "WSHPQFEK"             },  //  13
   { 7, "GSSGSSG"              },  //  14
   { 8, "TSLYKKAG"             },  //  15
   { 6, "HIEGRH"               },  //  16
   {11, "EQKLISEEDLN"          },  //  17
   { 8, "HHHHHGGS"             },  //  18
   { 5, "PPPPP"                },  //  19
   { 6, "EDEREE"               },  //  20
   {11, "AAALEHHHHHH"          },  //  21
   { 8, "GKTNFFEK"             },  //  22
   { 7, "GSRHHHH"              },  //  23
   {17, "HHHHHHSSGLEVLFQGP"    },  //  24
   { 5, "GGGGG"                },  //  25
   { 5, "SSTSS"                },  //  26
   { 5, "PPAPP"                },  //  27
   { 5, "APAPA"                },  //  28
   { 9, "ENLYFQGHM"            },  //  29
   { 5, "IDPFT"                },  //  30
   { 5, "KKKKK"                },  //  31
   { 5, "KSCDK"                },  //  32
   { 7, "EGKPIPN"              },  //  33
   { 5, "SSSVD"                },  //  34
   { 5, "MGRGS"                },  //  35
   { 8, "DCGCKPCI"             },  //  36
   { 7, "APAATGA"              },  //  37
   { 5, "KKKAA"                },  //  38
   { 5, "KKTSS"                },  //  39
   {11, "GGSGGGGSGGG"          },  //  40
   { 7, "HHHHSSG"              },  //  41
   {12, "PTTENLYFQGAM"         },  //  42
   { 6, "SSSSTQ"               },  //  43
   { 7, "EAQEEEE"              },  //  44
   { 7, "EEEEEEE"              },  //  45
   { 5, "AQAQE"                },  //  46
   { 6, "EPEEPA"               },  //  47
   { 7, "DDDDEDD"              },  //  48
   { 6, "QQQQQG"               },  //  49
   { 5, "PPAAP"                },  //  50
   { 7, "VKPEVKP"              },  //  51
   { 5, "PAATS"                },  //  52
   { 5, "SHMAS"                },  //  53
   { 6, "RRGKKK"               },  //  54
   { 7, "RGEGGFG"              },  //  55
   {10, "HHHHHHSQDP"           },  //  56
   { 5, "GPSSG"                },  //  57
   { 6, "GGAKRH"               },  //  58
   { 5, "APEDP"                },  //  59
   { 5, "PSPPP"                },  //  60
   { 6, "PSVSPS"               },  //  61
   { 7, "HHHHHMA"              },  //  62
   { 6, "DSVISS"               },  //  63
   { 6, "EPTESS"               },  //  64
   { 7, "AAVGGAA"              },  //  65
   { 5, "SPSPG"                },  //  66
   { 4, "GHMA"                 },  //  67
   { 5, "RSEED"                },  //  68
   { 7, "FPSPESS"              },  //  69
   { 5, "PAPPP"                },  //  70
   { 7, "QQQQQQQ"              },  //  71
   {11, "IKSHHNVGGLP"          },  //  72
   { 5, "EEEED"                },  //  73
   { 5, "YKDDD"                },  //  74
   { 6, "LDNGED"               },  //  75
   { 6, "PPPRPK"               },  //  76
   { 6, "GHHHHH"               },  //  77
   { 7, "NLREDGE"              },  //  78
   { 6, "PSPPSP"               },  //  79
   { 5, "YFQSM"                },  //  80
   {15, "GGLNDIFEAQKIEWH"      },  //  81
   { 6, "NGDTPS"               },  //  82
   { 7, "GRGHHHH"              },  //  83
   { 5, "GAMDP"                },  //  84
   { 5, "EAEED"                },  //  85
   { 7, "EDDEDED"              },  //  86
   { 5, "DDDKD"                },  //  87
   { 6, "QQREEG"               },  //  88
   { 6, "DIPESQ"               },  //  89
   { 5, "KKGKS"                },  //  90
   { 5, "RGEET"                },  //  91
   { 7, "GTHHHHH"              },  //  92
   { 9, "AAHHHHHHH"            },  //  93
   { 5, "ESSSS"                },  //  94
   { 5, "DAPDI"                },  //  95
   { 6, "ASIGQA"               },  //  96
   { 6, "TTTATT"               },  //  97
   { 5, "NNNNN"                },  //  98
   { 5, "RRRGR"                },  //  99
   { 6, "GQRKRR"               },  // 100
   { 5, "SMAEG"                },  // 101
   { 5, "SKKKK"                },  // 102
   { 6, "RPQLDS"               },  // 103
   { 6, "TSAETP"               },  // 104
   { 6, "RGRPRG"               },  // 105
   { 6, "DHSPAP"               },  // 106
   { 5, "AAPPA"                },  // 107
   { 7, "GGGSSSG"              },  // 108
   { 6, "QQQQQP"               },  // 109
   { 7, "AGAVAGG"              },  // 110
   { 6, "SMAAGG"               },  // 111
   { 5, "KKSKK"                },  // 112
   { 5, "REEEE"                },  // 113
   { 6, "SEFGSS"               },  // 114
   { 5, "EKKTE"                },  // 115
   { 5, "EKKNS"                },  // 116
   { 5, "GGGDD"                },  // 117
   { 8, "GEKHHHHH"             },  // 118
   { 6, "PPAPAG"               },  // 119
   { 5, "DLVPR"                },  // 120
   { 6, "RGSMAS"               },  // 121
   { 6, "GSETMA"               },  // 122
   { 7, "EEEKKKE"              },  // 123
   { 7, "RGGGGSG"              },  // 124
   { 6, "KKPKNK"               },  // 125
   { 6, "RSVRSN"               },  // 126
   { 5, "DEEDE"                },  // 127
   { 6, "AGEGPA"               },  // 128
   { 7, "SHHHHHH"              },  // 129
   { 6, "EDDESD"               },  // 130
   { 5, "NSSSS"                },  // 131
   { 5, "CGYSD"                },  // 132
   { 7, "SGSGGGS"              },  // 133
   { 5, "MEEEE"                },  // 134
   { 5, "GLVPR"                },  // 135
   { 5, "EKKKS"                },  // 136
   { 5, "GVPRG"                },  // 137
   { 6, "TDNGNS"               },  // 138
   { 7, "GMDELYK"              },  // 139
   { 6, "DEGHHH"               },  // 140
   { 5, "GGSRS"                },  // 141
   { 7, "GAHHHHH"              },  // 142
   { 7, "AMADIGS"              },  // 143
   { 7, "HHLHHHG"              },  // 144
   { 6, "GGKKKK"               },  // 145
   { 7, "SDEEDSS"              },  // 146
   { 5, "EEEEG"                },  // 147
   { 6, "SGDDDD"               },  // 148
   { 6, "AQSTSA"               },  // 149
   { 5, "PPPPQ"                },  // 150
   { 5, "GSMTD"                },  // 151
   {10, "GGGGSGGGGS"           },  // 152
   { 5, "EEDDD"                },  // 153
   { 5, "KKEKK"                },  // 154
   { 5, "STTST"                },  // 155
   { 8, "AELAAATA"             },  // 156
   { 6, "VDHHHH"               },  // 157
   { 8, "STSHHHHH"             },  // 158
   {21, "HHHHHHHHHSSGHIDDDDKHM"},  // 159
   { 6, "SDGKDD"               },  // 160
   { 8, "GSHMLEDP"             },  // 161
   { 6, "KSGYKD"               },  // 162
   { 5, "DEDSD"                },  // 163
   { 5, "DSDEE"                },  // 164
   { 6, "GGHNSS"               },  // 165
   { 5, "KSASS"                },  // 166
   { 5, "GSHGM"                },  // 167
   { 5, "MASPA"                },  // 168
   { 8, "SAWSHPQF"             },  // 169
   { 5, "GSEED"                },  // 170
   { 8, "ENLYFQGS"             }}; // 171



struct SP
  {
  double ac_elong[21];
  double a_init;
  double fs_n, fs_c;
  };

void SetPoten(int pat, int k)
{
static SP sp[3][K_MAX]={
  //     CYS,     MET,     PHE,     ILE,     LEU,     VAL,     TRP,     TYR,     ALA,     GLY,     THR,     SER,     GLN,     ASN,     GLU,     ASP,     HIS,    ARG,      LYS,    PRO,      UNK    a_init     fs_n     fs_c
 {{{  1.8445, -0.1274,  1.4886,  1.6735,  0.9861,  1.2418,  2.7407,  1.4576,  0.1118, -0.3110, -0.0091, -0.9648, -0.6926, -0.4590, -0.8077, -0.5130, -0.6628, -0.3163, -0.4954, -0.8322,  0.0000},  5.4604, -6.4977, -6.4968},   // 0 - Оптимальное Sw-0.02/t^2
  {{  1.2544, -0.0866,  1.0124,  1.1381,  0.6706,  0.8445,  1.8639,  0.9913,  0.0760, -0.2115, -0.0062, -0.6562, -0.4710, -0.3122, -0.5493, -0.3489, -0.4508, -0.2151, -0.3369, -0.5660,  0.0000},  3.7136, -4.4191, -4.4185},   // 1 - Оптимальное Sw50 && auc получено оптимизацией температуры
  {{  2.4188, -0.2387,  1.6671,  2.0799,  1.2732,  1.4630,  3.1954,  1.5785, -0.0639, -0.5250, -0.1529, -1.1092, -0.7109, -0.7780, -0.9353, -0.8225, -1.1264, -0.3427, -0.5974, -1.1507,  0.0000},  5.2510,  2.5911,  0.9429},   // 2 - Оценки усреднённые по позициям. Оптимальное Sw-0.02/t^2
  {{  1.6725, -0.1650,  1.1527,  1.4381,  0.8803,  1.0116,  2.2094,  1.0914, -0.0442, -0.3630, -0.1057, -0.7670, -0.4915, -0.5379, -0.6467, -0.5687, -0.7788, -0.2370, -0.4131, -0.7956,  0.0000},  3.6308,  1.7916,  0.6520},   // 3 - Оценки усреднённые по позициям. Оптимальное Sw50 && auc получено оптимизацией температуры
  {{  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000},  0,       0,       0,    }},  // 4 - Нулевые потенциалы

 {{{  1.7855, -0.1737,  1.4041,  1.6278,  0.9456,  1.2003,  2.6406,  1.4318,  0.1009, -0.3203, -0.0192, -0.9396, -0.6972, -0.4438, -0.7666, -0.5132, -0.0362, -0.3542, -0.4853, -0.8504,  0.0000},  5.4041, -6.4273, -6.4466},   // 0 - Оптимальное Sw-0.02/t^2
  {{  1.1684, -0.1137,  0.9188,  1.0652,  0.6188,  0.7854,  1.7279,  0.9369,  0.0660, -0.2096, -0.0126, -0.6148, -0.4562, -0.2904, -0.5016, -0.3358, -0.0237, -0.2318, -0.3176, -0.5565,  0.0000},  3.5363, -4.2058, -4.2185},   // 1 - Оптимальное Sw50 && auc получено оптимизацией температуры
  {{  2.3572, -0.1967,  1.5534,  1.9840,  1.1677,  1.3736,  3.0945,  1.5248, -0.0741, -0.5569, -0.1807, -1.0801, -0.6953, -0.7213, -0.9058, -0.8194, -0.0933, -0.3896, -0.5842, -1.1795,  0.0000},  5.3091,  3.0649,  1.1571},   // 2 - Оценки усреднённые по позициям. Оптимальное Sw-0.02/t^2
  {{  1.5956, -0.1331,  1.0515,  1.3430,  0.7904,  0.9298,  2.0947,  1.0322, -0.0502, -0.3770, -0.1223, -0.7311, -0.4707, -0.4883, -0.6132, -0.5547, -0.0632, -0.2637, -0.3955, -0.7984,  0.0000},  3.5938,  2.0747,  0.7833},   // 3 - Оценки усреднённые по позициям. Оптимальное Sw50 && auc получено оптимизацией температуры
  {{  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000},  0,       0,       0,    }},  // 4 - Нулевые потенциалы

 {{{  1.6559,  0.1355,  1.5356,  1.4983,  1.0375,  1.1006,  2.3728,  1.5854,  0.0436, -0.1797, -0.0344, -0.9096, -0.6062, -0.3939, -0.7278, -0.5486, -0.0166, -0.4009, -0.5519, -0.8824,  0.0000},  5.2553, -6.5076, -6.3774},   // 0 - Оптимальное Sw-0.02/t^2
  {{  1.0980,  0.0899,  1.0183,  0.9935,  0.6880,  0.7298,  1.5734,  1.0513,  0.0289, -0.1192, -0.0228, -0.6032, -0.4020, -0.2612, -0.4826, -0.3638, -0.0110, -0.2658, -0.3660, -0.5851,  0.0000},  3.4848, -4.3152, -4.2289},   // 1 - Оптимальное Sw50 && auc получено оптимизацией температуры
  {{  2.1239,  0.0697,  1.7384,  1.8623,  1.2839,  1.2617,  2.8732,  1.7961, -0.0715, -0.4354, -0.1518, -1.0900, -0.6191, -0.5098, -0.8548, -0.8930, -0.0795, -0.5002, -0.6018, -1.1907,  0.0000},  5.0979,  2.7147,  0.8467},   // 2 - Оценки усреднённые по позициям. Оптимальное Sw-0.02/t^2
  {{  1.3666,  0.0448,  1.1186,  1.1983,  0.8261,  0.8118,  1.8487,  1.1557, -0.0460, -0.2802, -0.0977, -0.7013, -0.3984, -0.3280, -0.5500, -0.5746, -0.0512, -0.3218, -0.3872, -0.7661,  0.0000},  3.2802,  1.7467,  0.5448},   // 3 - Оценки усреднённые по позициям. Оптимальное Sw50 && auc получено оптимизацией температуры
  {{  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000},  0,       0,       0,    }}}; // 4 - Нулевые потенциалы

 // 0 - Оптимальное Sw-0.02/t^2
 // 1 - Оптимальное Sw50 && auc получено оптимизацией температуры
 // 2 - Оценки усреднённые по позициям. Оптимальное Sw-0.02/t^2
 // 3 - Оценки усреднённые по позициям. Оптимальное Sw50 && auc получено оптимизацией температуры
 // 4 - Нулевые потенциалы

int i;
if(pat>2 || pat<0) pat=2;
PATTERN_T = pat;
if(pat==0) gpt=0;
if(pat==1) gpt=1;
if(pat==2) gpt=GPT;
if(k<0 || k>=K_MAX) k=K_MAX-1;
for(i=0; i<=20; i++) ac_elong[i] = sp[pat][k].ac_elong[i];
a_init = sp[pat][k].a_init;
fs_n   = sp[pat][k].fs_n;
fs_c   = sp[pat][k].fs_c;
}

double Unt_Pot(void)
{
int i,j;
double tunt;
tunt=0;
for(j=0; j<=20; j++) tunt += wac_elong[j] * ac_elong[j] * ac_elong[j];
tunt += wa_init * a_init * a_init;
tunt += wa_fs   * fs_n   * fs_n;
tunt += wa_fs   * fs_c   * fs_c;
tunt = sqrt(tunt);
return(tunt);
}

int Num_Ac(unsigned char a)
{
static const int ind[0x100] = {
  20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
  20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
  //   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V   W   X   Y   Z
  20,  8, 20,  0, 15, 14,  2,  9, 16,  3, 20, 18,  4,  1, 13, 20, 19, 12, 17, 11, 10, 20,  5,  6, 20,  7, 20, 20, 20, 20, 20, 20,
  20,  8, 20,  0, 15, 14,  2,  9, 16,  3, 20, 18,  4,  1, 13, 20, 19, 12, 17, 11, 10, 20,  5,  6, 20,  7, 20, 20, 20, 20, 20, 20,
  20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
  20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
  20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
  20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20};
return(ind[a]);
}


void SEQ_AC::Init(int n)
{
int i;
if(ac!=NULL) delete [] ac;
ac = new ACIDP[n];
num=0;
size=n;
for(i=0; i<n; i++) ac[i].reinit();
}

void SEQ_AC::ReInit(void)
{
int i;
for(i=0; i<num; i++) ac[i].reinit();
num=0;
}

void SEQ_AC::LoadFasta(char *seq_fasta)
{
int i, i0, n;
i0=0;
calc_pattern=0;
strcpy(info, "(no_name)");
if(*seq_fasta=='>')
  {
  for(i=0; seq_fasta[i]!=0; i++)
    {
    if(seq_fasta[i]=='\n')
      {
      i0=i+1;
      break;
      }
    }
  for(i=1; i<i0; i++)
    {
    if(seq_fasta[i]<=0x20 || (i+1)>=100) break;
    info[i-1] = seq_fasta[i];
    info[i] = 0;
    }
  }
for(i=i0,n=0; seq_fasta[i]!=0; i++)
  {
  if(seq_fasta[i]>0x20) n++;
  }
if(n>size) this->Init(n);
num = n;
for(i=i0,n=0; seq_fasta[i]!=0; i++)
  {
  if(seq_fasta[i]<=0x20) continue;
  ac[n].reinit();
  ac[n].ac = seq_fasta[i];
  ac[n].ind_ac = Num_Ac(seq_fasta[i]);
  n++;
  }
}

void SEQ_AC::LoadFastaFile(const char *name)
{
FILE *inp;
char *s, ss;
int n;
calc_pattern=0;
inp = fopen(name, "r");
if(inp == NULL)
  {
  fprintf(stderr, "Cannot open \"%s\"\n", name);
  exit(errno);
  }
for(n=0,ss=fgetc(inp); !feof(inp); ss=fgetc(inp)) n++;
s = new char[n+2];
rewind(inp);
for(n=0,ss=fgetc(inp); !feof(inp); ss=fgetc(inp))
  {
  s[n] = ss;
  n++;
  s[n]=0;
  }
this->LoadFasta(s);
delete [] s;
fclose(inp);
}

double E_Add(double a, double b)
{
double rz;
if(b<a) return(E_Add(b,a));
b-=a;
if(b>100) return(a);
rz = a - log(1. + exp(-b));
return(rz);
}

void Check_Pattern(SEQ_AC &sa, int i0, int im, PATTERN &pt)
{
int i,j, k;
if(im-i0 < pt.n) return; // i=i0;i<im
j=0;      i=i0+j; if(sa[i].ac != pt[j]) return;
j++;      i=i0+j; if(sa[i].ac != pt[j]) return;
j=pt.n-1; i=i0+j; if(sa[i].ac != pt[j]) return;
j--;      i=i0+j; if(sa[i].ac != pt[j]) return;
for(k=0,i=i0,j=0; j<pt.n; i++,j++)
  {
  if(pt[j] == 'X') continue;
  if(sa[i].ac == pt[j])
       k++;
  else k-=CH_PATTERN;
  if(k<-pt.n) return;
  }
//if(k<pt.n) return;
if(k<=0) return;
//for(i=i0,j=0; j<pt.n; i++,j++) sa[i].ept=E_PATTERN;
im = i0+pt.n;
//if(im <= 40) i0=0;
//if(i0+40 >= Length(sa)) im=Length(sa);
if(i0 < 40) i0=0;
if(Length(sa)-im < 40) im=Length(sa);
for(i=i0; i<im; i++) sa[i].ept=E_PATTERN;
}


void Pred_Pattern(SEQ_AC &sa)
{
int i0, im, l, i;
for(i=0; i<Length(sa); i++) sa[i].ept=0;
sa.calc_pattern = 1;
if(PATTERN_T==0)
  {
  sa.calc_pattern = 1;
  return;
  }
im = Length(sa);
for(l=0; l<gpt; l++) for(i=0; i<im; i++)
  {
  Check_Pattern(sa, i, im, pt_all[l]);
  }
}

void PredictNC(SEQ_AC &sa)
{
int i, ia;
if(Length(sa)<1) return;
i = 0;
ia = sa[i].ind_ac;
sa[i].e3d_n = sa[i].ept;
sa[i].elp_n = unt * (ac_elong[ia] + fs_n);
for(i=1; i<Length(sa); i++)
  {
  ia  = sa[i  ].ind_ac;
  sa[i].e3d_n = E_Add(sa[i-1].e3d_n,  sa[i-1].elp_n + unt * a_init) + sa[i].ept;
  sa[i].elp_n = E_Add(sa[i-1].elp_n,  sa[i-1].e3d_n + unt * a_init) + unt * ac_elong[ia];
  }
}

void PredictCN(SEQ_AC &sa)
{
int i, iap;
if(Length(sa)<1) return;
i = Length(sa)-1;
sa[i].e3d_c = 0;
sa[i].elp_c = unt * fs_c;
for(i=Length(sa)-2; i>=0; i--)
  {
  iap = sa[i+1].ind_ac;

  sa[i].e3d_c = E_Add(sa[i+1].e3d_c + sa[i+1].ept, sa[i+1].elp_c + unt*(a_init + ac_elong[iap]));
  sa[i].elp_c = E_Add(sa[i+1].elp_c + unt * ac_elong[iap],    sa[i+1].e3d_c + unt*a_init + sa[i+1].ept);
  }
}

void ACIDP::CalcP(void)
{
double e3d, elp, sp;
e3d = e3d_n + e3d_c;
elp = elp_n + elp_c;
p3d = 1.;
plp = exp(-elp+e3d);
sp = p3d+plp;
p3d/=sp;
plp/=sp;
}

void Predict(SEQ_AC &sa)
{
int i;
//printf("Predict %s %d\n", sa.info, sa.calc_pattern);
if(sa.calc_pattern==0) Pred_Pattern(sa);
PredictNC(sa);
PredictCN(sa);
for(i=0; i<Length(sa); i++) sa[i].CalcP();
}

void Predict_sh(SEQ_AC &sa)
{
int i;
Pred_Pattern(sa);
for(i=0; i<Length(sa); i++)
  {
  if(sa[i].ept<0.5)
    {
    sa[i].p3d = 0.5;
    sa[i].plp = 0.5;
    }
  else
    {
    sa[i].p3d = 0;
    sa[i].plp = 1;
    }
  }
}


