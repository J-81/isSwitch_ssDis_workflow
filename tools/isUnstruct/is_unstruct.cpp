#include "pred_del.hpp"
#include <ctype.h>

const char *release = "2.02";
char file_name[0x200]="";
int use_pattern=2,
    opt_type=1;      // for autors of IsUnstruct only
char use_pattern_info[3][40]= {
  "Patterns not used",
  "Only H6 pattern was used",
  "All patterns were used"};
int long_disp=0, short_disp=0, file_save=1;

void CreateOutName(char *out_name, const char *inp_name, const char *ext)
{
int i0, i1, i;
i1=i0=0;
for(i=0; inp_name[i]!=0; i++)
  {
  if(inp_name[i]=='.') i1=i;
  if(inp_name[i]=='/' || inp_name[i]=='\\' )
    {
    i0=i+1;
    i1=i+1;
    }
  }
if(i1<=i0) i1=strlen(inp_name);
strncpy(out_name, inp_name+i0, i1-i0); out_name[i1-i0]=0;
if(*ext!='.') strcat(out_name, ".");
strcat(out_name, ext);
//printf("\"%s\" < \"%s\" + \"%s\"\n", out_name, inp_name, ext);
}

void SaveLong(SEQ_AC &sa, char *inp_name)
{
char *out_name;
int i;
FILE *out;
char c, pt;
out_name = new char[strlen(inp_name)+5];
CreateOutName(out_name, inp_name, "iul");
out = fopen(out_name, "w");
if(out==NULL)
  {
  fprintf(stderr, "Cannot create %s\n", out_name);
  return;
  }
fprintf(out, "# IsUnstruct %s release / long format\n", release);
fprintf(out, "# %s\n", use_pattern_info[use_pattern]);
if(opt_type!=1) fprintf(out, "# opt_type=%d. The prediction can be incorrect!\n", opt_type);
fprintf(out, "# %s\n\n", sa.info);
for(i=0; i<Length(sa); i++)
  {
  pt=' ';
  if(sa[i].ept > (E_PATTERN/2)) pt='P';
  if(sa[i].plp>0.5) c='U';
  else              c='s';
  fprintf(out, "%4d %c %c %c %5.3f\n", i+1, sa[i].ac, c, pt, sa[i].plp);
  }
fclose(out);
delete [] out_name;
}

void SaveShort(SEQ_AC &sa, char *inp_name)
{
char *out_name;
int i, j, i0, np;
FILE *out;
char c, q, cpt, tp, seq[120], qp[120], prd[120], pt[120];
out_name = new char[strlen(inp_name)+5];
CreateOutName(out_name, inp_name, "ius");
out = fopen(out_name, "w");
if(out==NULL)
  {
  fprintf(stderr, "Cannot create %s\n", out_name);
  return;
  }
fprintf(out, "# IsUnstruct %s release / short format\n", release);
fprintf(out, "# %s\n", use_pattern_info[use_pattern]);
if(opt_type!=1) fprintf(out, "# opt_type=%d. The prediction can be incorrect!\n", opt_type);
fprintf(out, "# %s\n\n", sa.info);
*seq=*qp=*prd=*pt=0;
for(i0=0; i0<Length(sa); i0+=100)
  {
  np=0;
  for(i=i0, j=0; i<(i0+100); i++)
    {
    if(i < Length(sa))
      {
      cpt = ' ';
      if(sa[i].ept > (E_PATTERN/2))
        {
        np++;
        cpt='P';
        }
      c = sa[i].ac;
      if(sa[i].plp>0.5) tp='U';
      else              tp='s';
      q='0';
      if(sa[i].plp>=0.1  ) q='1';
      if(sa[i].plp>=0.2  ) q='2';
      if(sa[i].plp>=0.3  ) q='3';
      if(sa[i].plp>=0.4  ) q='4';
      if(sa[i].plp>=0.5  ) q='5';
      if(sa[i].plp>=0.6  ) q='6';
      if(sa[i].plp>=0.7  ) q='7';
      if(sa[i].plp>=0.8  ) q='8';
      if(sa[i].plp>=0.9  ) q='9';
      //if(sa[i].plp>=0.995) q='U';
      }
    else
      {
      cpt= ' ';
      tp = ' ';
      c  = ' ';
      q  = ' ';
      }
    if((i%100)!=0 && (i%10)==0)
      {
      pt[j] =' ';
      prd[j]=' ';
      seq[j]=' ';
      qp [j]=' ';
      j++;
      }
    pt[j]  =cpt;   pt[j+1]  = 0;
    prd[j] = tp;   prd[j+1] = 0;
    seq[j] =  c;   seq[j+1] = 0;
    qp[j]  =  q;   qp[ j+1] = 0;
    j++;
    }
  if(i0+100 <= Length(sa))
    {
    fprintf(out, "sequence   %5d %-109s %4d\n",   i0+1, seq, i0+100);
    fprintf(out, "state      %5d %-109s %4d\n",   i0+1, prd, i0+100);
    if(np>0) fprintf(out, "pattern    %5d %-109s %4d\n", i0+1, pt,  i0+100);
    fprintf(out, "probability%5d %-109s %4d\n\n", i0+1, qp,  i0+100);
    }
  else
    {
    fprintf(out, "sequence   %5d %-109s %4d\n",   i0+1, seq, Length(sa));
    fprintf(out, "state      %5d %-109s %4d\n",   i0+1, prd, Length(sa));
    if(np>0) fprintf(out, "pattern    %5d %-109s %4d\n", i0+1, pt, Length(sa));
    fprintf(out, "probability%5d %-109s %4d\n\n", i0+1, qp,  Length(sa));
    }
  }
fclose(out);
delete [] out_name;
}

void SaveLong(SEQ_AC &sa)
{
int i;
char c, pt;
fprintf(stdout, "# IsUnstruct %s release / long format\n", release);
fprintf(stdout, "# %s\n", use_pattern_info[use_pattern]);
if(opt_type!=1) fprintf(stdout, "# opt_type=%d. The prediction can be incorrect!\n", opt_type);
fprintf(stdout, "# %s\n\n", sa.info);
for(i=0; i<Length(sa); i++)
  {
  pt=' ';
  if(sa[i].ept > (E_PATTERN/2)) pt='P';
  if(sa[i].plp>0.5) c='U';
  else              c='s';
  fprintf(stdout, "%4d %c %c %c %5.3f\n", i+1, sa[i].ac, c, pt, sa[i].plp);
  }
}

void SaveShort(SEQ_AC &sa)
{
int i, j, i0, np;
char c, q, cpt, tp, seq[120], qp[120], prd[120], pt[120];
fprintf(stdout, "# IsUnstruct %s release / short format\n", release);
fprintf(stdout, "# %s\n", use_pattern_info[use_pattern]);
if(opt_type!=1) fprintf(stdout, "# opt_type=%d. The prediction can be incorrect!\n", opt_type);
fprintf(stdout, "# %s\n\n", sa.info);
*seq=*qp=*prd=*pt=0;
for(i0=0; i0<Length(sa); i0+=100)
  {
  np=0;
  for(i=i0, j=0; i<(i0+100); i++)
    {
    if(i < Length(sa))
      {
      cpt = ' ';
      if(sa[i].ept > (E_PATTERN/2))
        {
        np++;
        cpt='P';
        }
      c = sa[i].ac;
      if(sa[i].plp>0.5) tp='U';
      else              tp='s';
      q='0';
      if(sa[i].plp>=0.1  ) q='1';
      if(sa[i].plp>=0.2  ) q='2';
      if(sa[i].plp>=0.3  ) q='3';
      if(sa[i].plp>=0.4  ) q='4';
      if(sa[i].plp>=0.5  ) q='5';
      if(sa[i].plp>=0.6  ) q='6';
      if(sa[i].plp>=0.7  ) q='7';
      if(sa[i].plp>=0.8  ) q='8';
      if(sa[i].plp>=0.9  ) q='9';
      //if(sa[i].plp>=0.995) q='U';
      }
    else
      {
      cpt= ' ';
      tp = ' ';
      c  = ' ';
      q  = ' ';
      }
    if((i%100)!=0 && (i%10)==0)
      {
      pt[j] =' ';
      prd[j]=' ';
      seq[j]=' ';
      qp [j]=' ';
      j++;
      }
    pt[j]  =cpt;   pt[j+1]  = 0;
    prd[j] = tp;   prd[j+1] = 0;
    seq[j] =  c;   seq[j+1] = 0;
    qp[j]  =  q;   qp[ j+1] = 0;
    j++;
    }
  if(i0+100 <= Length(sa))
    {
    fprintf(stdout, "sequence   %5d %-109s %4d\n",   i0+1, seq, i0+100);
    fprintf(stdout, "state      %5d %-109s %4d\n",   i0+1, prd, i0+100);
    if(np>0) fprintf(stdout, "pattern    %5d %-109s %4d\n", i0+1, pt,  i0+100);
    fprintf(stdout, "probability%5d %-109s %4d\n\n", i0+1, qp,  i0+100);
    }
  else
    {
    fprintf(stdout, "sequence   %5d %-109s %4d\n",   i0+1, seq, Length(sa));
    fprintf(stdout, "state      %5d %-109s %4d\n",   i0+1, prd, Length(sa));
    if(np>0) fprintf(stdout, "pattern    %5d %-109s %4d\n", i0+1, pt, Length(sa));
    fprintf(stdout, "probability%5d %-109s %4d\n\n", i0+1, qp,  Length(sa));
    }
  }
}

void LoadKey(int argc, char **argv)
{
int j;
for(j=1; j<argc; j+=2)
  {
  if(argv[j][0] == '-')
    {
    switch(argv[j][1])
      {
      case 'u': use_pattern = atoi(argv[j+1]); break;
      case 'o': opt_type    = atoi(argv[j+1]); break; // for autors of IsUnstruct only
      case 'l': long_disp   = atoi(argv[j+1]); break;
      case 's': short_disp  = atoi(argv[j+1]); break;
      case 'f': file_save   = atoi(argv[j+1]); break;
      }
    }
  else
    {
    strcpy(file_name, argv[j]);
    break;
    }
  }
SetPoten(use_pattern, opt_type);
}

int main(int argc, char **argv)
{
SEQ_AC sa;
LoadKey(argc, argv);          // printf("LoadKey %s\n", file_name);
if(file_name[0]==0)
  {
  printf("Usage: ./IsUnstruct [-use_pattern 0,1,2] [-long_disp 0,1] [-short_disp 0,1] [-file_save 0,1] <file_in_fasta_format>\n");
  printf("-use_pattern -- 0           Without patterns\n");
  printf("                1           With HHHHHH pattern only\n");
  printf("                2 (default) With all patterns\n");
  printf("-long_disp   -- 0 (default) Do not show long  output by the terminal\n");
  printf("             -- 1                  Show long  output by the terminal\n");
  printf("-short_disp  -- 0 (default) Do not show short output by the terminal\n");
  printf("             -- 1                  Show short output by the terminal\n");
  printf("-file_save   -- 0           Do not save output files\n");
  printf("             -- 1 (default)        Save output files\n");
  return(0);
  }
sa.LoadFastaFile(file_name);  // printf("Load %s\n", file_name);
Predict(sa);                  // printf("Predict\n");
if(long_disp >0) SaveLong( sa);
if(short_disp>0) SaveShort(sa);
if(file_save>0)
  {
  SaveLong( sa, file_name);     // printf("SaveLong\n");
  SaveShort(sa, file_name);     // printf("SaveShort\n");
  }
return(0);
}
