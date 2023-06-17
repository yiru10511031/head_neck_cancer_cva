LIBNAME TMP1 'D:\Downloads';
/*Cor_HNCA_RADI_CHEMO_CVA���w�X�֩�����A�L�{���b�o��*/
***Step6: ��ƾ�z&�έp���R;
*6-1: �إ�hn�ܼơA�N�e�f�����Y�V���P��|��;
DATA TMP1.CVAALL;
SET TMP1.Cor_HNCA_RADI_CHEMO_CVA;
attrib HN length=$10;
IF ACODE_ICD9_1 IN:('147' ) OR 
    ACODE_ICD9_2 IN:('147')  OR
	ACODE_ICD9_3 IN:('147' ) 
THEN HN='NPC'; 
ELSE IF ACODE_ICD9_1 IN:('140' '141' '143' '144' '145' '146' '148' '149') OR 
    ACODE_ICD9_2 IN:('140' '141' '143' '144' '145' '146'  '148' '149') OR
	ACODE_ICD9_3 IN:('140' '141' '143' '144' '145' '146'  '148' '149')
THEN HN='HNC'; 
ELSE HN=.;
RUN;
*8985;
PROC FREQ DATA=TMP1.CVAALL;
TABLE ID_SEX*HN;
RUN;

*6-2: �]���ʧO��U���ҥ~�A�٦���|�ȡA����J�Q��;
DATA TMP1.CVAALL;
SET TMP1.CVAALL;
IF ID_SEX='U' THEN DELETE;
IF ID_SEX='' THEN DELETE;
IF RADI='.' THEN RADI=0;
IF CHEMO='.' THEN CHEMO=0;
RUN;
*8977;

*6-3: �p��~��;
DATA TMP1.CVAALL;
SET TMP1.CVAALL;
birth_year= SUBSTR(ID_BIRTHDAY,1,4); *�]���榡�ݰ_�ӥu���~��S����A�B�D����榡�A���e�|�X���X�ͦ~;
MED_DAY=INPUT(FUNC_DATE, YYMMDD8.);
FORMAT MED_DAY YYMMDD10.;
AGE=round((MED_DAY-birth_year)/365.25, 0.1); **��X�N�媺������X���A�|�ˤ��J��p�ƲĤ@��;
RUN;

*6-4: �s���ܶ�;
*6-4-1: �~��*CVA;
PROC MEANS DATA=TMP1.CVAALL Maxdec=2 n mean std var min q1 median q3 max ;
CLASS hn CVA;
VAR AGE;
RUN;
proc sort data=TMP1.CVA_ALL;
by age;
run;
*�e��: ��|���B�Y�V���@�_;
*�����;
proc sgpanel data=TMP1.CVAALL ;
panelby hn / novarname sparse columns = 2 rows = 1; /*�T���C*/
*rowaxis label="Percentage (%)"; /*���O���wY, X�y�жb���W��*/
*colaxis label = "AGE";
  histogram age / group=cva transparency=0.5 binwidth=1;       /* SAS 9.4m2 */
  density age / type=kernel group=cva ; /* overlay density estimates */
run;
*Box-plot;
proc sgpanel data=TMP1.CVAALL ;
panelby hn / novarname sparse columns = 2 rows = 1; /*�T���C*/
*rowaxis label="Percentage (%)"; /*���O���wY, X�y�жb���W��*/
*colaxis label = "AGE";
hbox age / group=cva;
run;
*�~�֤���;
data TMP1.CVAALL ;
set TMP1.CVAALL ;
if age<40 then age_g='30-39';
else if age>=40 & age<50 then age_g='40-49';
else age_g=.;
run;

*6-5: ���s�p��CVA_YR��CVA_YR2;
DATA TMP1.CVAALL ;
SET TMP1.CVAALL ;
cva_yr2=intck('day', index, brain_index)/365.25;
RUN;

*6-6: ����;
*��|���Ӯ�;
DATA TMP1.CVA147;
SET TMP1.CVAALL ;
IF HN='NPC'; 
RUN;
*3944;
PROC FREQ DATA=TMP1.CVA147;
TABLES DM Liver MI Lip HTN  CVA;
RUN;

*��L�Y�V���Ӯ�;
DATA TMP1.CVA_ELSE;
SET TMP1.CVAALL ;
IF HN='HNC'; 
RUN;
*5033;
PROC FREQ DATA=TMP1.CVA_ELSE;
TABLES DM Liver MI Lip HTN  CVA;
RUN;


*6-7: ���O�ܶ�;
*6-7-1: �[��CVA��;
PROC FREQ DATA=TMP1.CVAALL;
TABLES hn*CVA;
RUN;

*6-7-2: �d��/�O��;
*1.��|��;
PROC FREQ DATA=TMP1.CVA147;
TABLES  CVA*ID_SEX  CVA*AGE_G CVA*dm  CVA*Liver 
CVA*MI  CVA*Lip  CVA*HTN  CVA*RADI   CVA*CHEMO/ expected chisq fisher;
RUN;
*2.�Y�V��;
PROC FREQ DATA=TMP1.CVA_ELSE;
TABLES  CVA*ID_SEX  CVA*AGE_G CVA*dm  CVA*Liver 
CVA*MI  CVA*Lip  CVA*HTN  CVA*RADI   CVA*CHEMO/ expected chisq or;
RUN;


*6-8: KM�s�����u Kaplan-Meier Method;
*1.��|��;
*�H�ʧO������;
proc lifetest data=TMP1.CVA147 plot=survival(test nocensor) ;
time cva_yr2*cva(0); 
STRATA id_SEX;
run;
*�H�~�֬�����;
proc lifetest data=TMP1.CVA147 plot=survival(test nocensor) ;
time cva_yr2*cva(0); 
STRATA age_g;
run;

*2.�Y�V��;
*�H�ʧO������;
proc lifetest data=TMP1.CVA_ELSE plot=survival(test nocensor) ;
time cva_yr2*cva(0); 
STRATA id_SEX;
run;
*�H�~�֬�����;
proc lifetest data=TMP1.CVA_ELSE plot=survival(test nocensor) ;
time cva_yr2*cva(0); 
STRATA age_g;
run;


*Cox proportional hazard regression model;
*1.��|��;
*cva_yr2����;
/*data TMP1.CVA147;
set TMP1.CVA147;
if cva_yr2>=2 & cva_yr2<4 then cvayr2_g='2<=yr<4';
else if cva_yr2<2 then cvayr2_g='<2';
else if cva_yr2>=4 & cva_yr2<6 then cvayr2_g='4<=yr<6';
else if cva_yr2>=6 & cva_yr2<8 then cvayr2_g='6<=yr<8';
else if cva_yr2>=8 then cvayr2_g='yr>=8';
else cvayr2_g=.;
run;*/


proc phreg data=TMP1.CVA147 plots(overlay)=survival;
class id_sex(param=ref ref='F') ;
model cva_yr2*cva(0)=id_sex/ RISKLIMITS;
run;
proc phreg data=TMP1.CVA147 plots(overlay)=survival;
class  age_g(param=ref ref='30-39') ;
model cva_yr2*cva(0)=age_g/ RISKLIMITS;
run;
/*proc phreg data=TMP1.CVA147 plots(overlay)=survival;
class cvayr2_g(param=ref ref='<2');
model cva_yr2*cva(0)=cvayr2_g/ RISKLIMITS;
run;*/
proc phreg data=TMP1.CVA147 plots(overlay)=survival;
class dm(ref='0');
model cva_yr2*cva(0)=dm/ RISKLIMITS;
run;
proc phreg data=TMP1.CVA147 plots(overlay)=survival;
class Liver(ref='0');
model cva_yr2*cva(0)= Liver/ RISKLIMITS;
run;
proc phreg data=TMP1.CVA147 plots(overlay)=survival;
class MI(ref='0');
model cva_yr2*cva(0)=MI / RISKLIMITS;
run;
proc phreg data=TMP1.CVA147 plots(overlay)=survival;
class Lip(ref='0');
model cva_yr2*cva(0)=Lip/ RISKLIMITS;
run;
proc phreg data=TMP1.CVA147 plots(overlay)=survival;
class HTN(ref='0') ;
model cva_yr2*cva(0)=HTN / RISKLIMITS;
run;
proc phreg data=TMP1.CVA147 plots(overlay)=survival;
class RADI(ref='0');
model cva_yr2*cva(0)= RADI/ RISKLIMITS;
run;
proc phreg data=TMP1.CVA147 plots(overlay)=survival;
class CHEMO(ref='0');
model cva_yr2*cva(0)= CHEMO/ RISKLIMITS;
run;

*2.�Y�V��;
*cva_yr2����;
/*data TMP1.CVA_ELSE;
set TMP1.CVA_ELSE;
if cva_yr2>=2 & cva_yr2<4 then cvayr2_g='2<=yr<4';
else if cva_yr2<2 then cvayr2_g='<2';
else if cva_yr2>=4 & cva_yr2<6 then cvayr2_g='4<=yr<6';
else if cva_yr2>=6 & cva_yr2<8 then cvayr2_g='6<=yr<8';
else if cva_yr2>=8 then cvayr2_g='yr>=8';
else cvayr2_g=.;
run;*/

proc phreg data=TMP1.CVA_ELSE plots(overlay)=survival;
class id_sex(param=ref ref='F') ;
model cva_yr2*cva(0)=id_sex/ RISKLIMITS;
run;
proc phreg data=TMP1.CVA_ELSE plots(overlay)=survival;
class  age_g(param=ref ref='30-39') ;
model cva_yr2*cva(0)=age_g/ RISKLIMITS;
run;

/*proc phreg data=TMP1.CVA_ELSE plots(overlay)=survival;
class cvayr2_g(param=ref ref='<2');
model cva_yr2*cva(0)=cvayr2_g/ RISKLIMITS;
run;*/
proc phreg data=TMP1.CVA_ELSE plots(overlay)=survival;
class dm(ref='0');
model cva_yr2*cva(0)=dm/ RISKLIMITS;
run;
proc phreg data=TMP1.CVA_ELSE plots(overlay)=survival;
class Liver(ref='0');
model cva_yr2*cva(0)= Liver/ RISKLIMITS;
run;
proc phreg data=TMP1.CVA_ELSE plots(overlay)=survival;
class MI(ref='0');
model cva_yr2*cva(0)=MI / RISKLIMITS;
run;
proc phreg data=TMP1.CVA_ELSE plots(overlay)=survival;
class Lip(ref='0');
model cva_yr2*cva(0)=Lip/ RISKLIMITS;
run;
proc phreg data=TMP1.CVA_ELSE plots(overlay)=survival;
class HTN(ref='0') ;
model cva_yr2*cva(0)=HTN / RISKLIMITS;
run;
proc phreg data=TMP1.CVA_ELSE plots(overlay)=survival;
class RADI(ref='0');
model cva_yr2*cva(0)= RADI/ RISKLIMITS;
run;
proc phreg data=TMP1.CVA_ELSE plots(overlay)=survival;
class CHEMO(ref='0');
model cva_yr2*cva(0)= CHEMO/ RISKLIMITS;
run;


proc phreg data=TMP1.CVA_ELSE plots(overlay)=survival;
class id_sex(param=ref ref='F') age_g(param=ref ref='30-39') cvayr2_g(param=ref ref='<2')
		dm(ref='0')  Liver(ref='0') MI(ref='0')  Lip(ref='0') HTN(ref='0') 
		RADI(ref='0') CHEMO(ref='0');
model cva_yr2*cva(0)=id_sex age_g cvayr2_g dm Liver MI Lip HTN RADI CHEMO/ RISKLIMITS;
run;
