#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#define _CRT_SECURE_NO_WARNINGS
#define pcm 10        //对于蒙卡程序，有自带的统计偏差，因此很难做到两个模型的k完全相等，
//所以设置3pcm作为偏差限，两个模型k偏差在3pcm以内，认为等效
#define inf 100000
#define cyc_max 20   //循环上限
#define NUM 3672    //颗粒数
//#define NUM 2902    //颗粒数
#define pi 3.1415926
#define h 4
#define avo 6.02E+23  //阿伏伽德罗常数
#define barn 1.0E+24

FILE* fp1;
FILE* fp2;
//char *RPT_filein;
//char *RPT_fileout;
double RPT_R;
double R1, R2, R3;
double deltK1, deltK2, deltK3;
double K_fcm;
double bias_fcm;
double K_rpt1, K_rpt2, K_rpt3;
double bias_rpt1, bias_rpt2, bias_rpt3;
double ZL_U5;
double ZL_U8;
double ZL_N14;
double ZL_Si;
double ZL_Si28;
double ZL_Si29;
double ZL_Si30;
double ZL_C;
double r_fuel;
double r_buf;
double r_pyc1;
double r_sic;
double r_pyc2;
double Rou_h2o;
double nd_u8;
double nd_u5;
double nd_n14;
double nd_buf;
double nd_pyc;
double nd_sic;
double nd_h;
double nd_o;
double L;
double r_fcm;
double r_gap;
double r_clad;

//读取输入卡信息
void Read_para()
{
	char none[500];
	FILE* fp;
	fp = fopen("Parameters", "r");

	fscanf(fp, "%s", none);

	fscanf(fp, "%s", none);
	fscanf(fp, "%lf", &r_fuel);//fuel    0.0380
	fscanf(fp, "%s", none);
	fscanf(fp, "%lf", &r_buf);//buffer  0.0400
	fscanf(fp, "%s", none);
	fscanf(fp, "%lf", &r_pyc1);//PyC1     0.0435
	fscanf(fp, "%s", none);
	fscanf(fp, "%lf", &r_sic);//SiC     0.0470
	fscanf(fp, "%s", none);
	fscanf(fp, "%lf", &r_pyc2);//PyC2     0.0490

	fscanf(fp, "%s", none);
	fscanf(fp, "%s", none);
	fscanf(fp, "%lf", &r_fcm);//fuel  0.600
	fscanf(fp, "%s", none);
	fscanf(fp, "%lf", &r_gap);//gap  0.606
	fscanf(fp, "%s", none);
	fscanf(fp, "%lf", &r_clad);//clad 0.691
	fscanf(fp, "%s", none);
	fscanf(fp, "%lf", &L);//water  0.824
	printf("L=%f\n", L);

	fscanf(fp, "%s", none);
	fscanf(fp, "%s", none);
	fscanf(fp, "%lf", &nd_u5);//92235.03c   2.2519E-03
	fscanf(fp, "%s", none);
	fscanf(fp, "%lf", &nd_u8);//92238.03c   3.1984E-02
	fscanf(fp, "%s", none);
	fscanf(fp, "%lf", &nd_n14);// 7014.03c   3.4236E-02
	fscanf(fp, "%s", none);
	fscanf(fp, "%lf", &nd_buf);//buffer  5.2675E-2
	fscanf(fp, "%s", none);
	fscanf(fp, "%lf", &nd_pyc);//pyc   9.5317E-2
	fscanf(fp, "%s", none);
	fscanf(fp, "%lf", &nd_sic);//SiC_SiC   4.7859e-2
	printf("nd_sic=%E\n", nd_sic);
	fscanf(fp, "%s", none);
	fscanf(fp, "%lf", &Rou_h2o);//Rou_h2o  1.0
	fclose(fp);
	nd_o = Rou_h2o * avo / 18.0 / barn;
	nd_h = 2 * nd_o;



}

//均匀化换算
void CalZL()
{
	/*r_fuel = 0.035;
	r_buf = 0.040;
	r_pyc1 = 0.0435;
	r_sic = 0.047;
	r_pyc2 = 0.049;*/


	ZL_U5 = NUM * 4.0 * pi * r_fuel * r_fuel * r_fuel * nd_u5 / 3.0;//计算U235的核子数
	ZL_U8 = NUM * 4.0 * pi * r_fuel * r_fuel * r_fuel * nd_u8 / 3.0;//计算U238的核子数
	ZL_N14 = NUM * 4.0 * pi * r_fuel * r_fuel * r_fuel * nd_n14 / 3.0;//计算N14的核子数
	ZL_Si = (NUM * 4.0 * pi * (r_sic * r_sic * r_sic - r_pyc1 * r_pyc1 * r_pyc1) / 3.0 + (pi * r_fcm * r_fcm * h - NUM * 4.0 * pi * r_pyc2 * r_pyc2 * r_pyc2 / 3.0)) * nd_sic;/*计算SiC的原子数由两部分组成
			颗粒内的SiC层             										基质里面抛开燃料颗粒的SiC													  */
	ZL_C = ZL_Si + (NUM * 4.0 * pi * (r_pyc1 * r_pyc1 * r_pyc1 - r_buf * r_buf * r_buf) / 3.0 + NUM * 4.0 * pi * (r_pyc2 * r_pyc2 * r_pyc2 - r_sic * r_sic * r_sic) / 3.0) * nd_pyc;
	ZL_C = ZL_C + (NUM * 4.0 * pi * (r_buf * r_buf * r_buf - r_fuel * r_fuel * r_fuel) / 3.0) * nd_buf;//计算所有的C

	printf("ZL_Si=%E\n", ZL_Si);
	printf("ZL_C=%E\n", ZL_C);
}
//颗粒模型输入卡自动生成
void FCM_CARD()
{
	FILE* fp;
	fp = fopen("FCM", "w");
	fprintf(fp, " \n");
	fprintf(fp, "set title \"HTGR Depletion Benchmark, Pebble Bed Model\"\n");
	fprintf(fp, "\n");
	fprintf(fp, "%% ---------------------------------------------------------------------------- - \n");
	fprintf(fp, "\n");
	fprintf(fp, "%% Geometry:\n");
	fprintf(fp, "\n");
	fprintf(fp, "%% ---------------------------------------------------------------------------- - \n");
	fprintf(fp, "\n");
	fprintf(fp, "%% -- - Coated particle : \n");
	fprintf(fp, "\n");
	fprintf(fp, "particle 1\n");
	fprintf(fp, "\n");
	fprintf(fp, "fuel    %.4f\n", r_fuel);
	fprintf(fp, "buffer  %.4f\n", r_buf);
	fprintf(fp, "PyC     %.4f\n", r_pyc1);
	fprintf(fp, "SiC     %.4f\n", r_sic);
	fprintf(fp, "PyC     %.4f\n", r_pyc2);
	fprintf(fp, "SiC\n");
	fprintf(fp, "\n");
	fprintf(fp, "%% -- - Universe surrounding particles : \n");
	fprintf(fp, "\n");
	fprintf(fp, "surf 10 inf\n");
	fprintf(fp, "\n");
	fprintf(fp, "cell 10 2 SiC -10\n");
	fprintf(fp, "\n");
	fprintf(fp, "%% -- - Read particles : \n");
	fprintf(fp, "\n");
	fprintf(fp, "pbed 10 2 \"2.inp\"\n");
	fprintf(fp, "\n");
	fprintf(fp, "%% -- - Pebble\n");
	fprintf(fp, "\n");
	fprintf(fp, "surf 1 cylz   0.0 0.0 %.4f -2 2\n", r_fcm);
	fprintf(fp, "surf 2 cylz   0.0 0.0 %.4f -2 2\n", r_gap);
	fprintf(fp, "surf 3 cylz   0.0 0.0 %.4f -2 2\n", r_clad);
	fprintf(fp, "surf 4 cuboid %.4f %.4f %.4f %.4f -2 2\n", -1.0 * L, L, -1.0 * L, L);
	fprintf(fp, "\n");
	fprintf(fp, "cell 1 0 fill 10 -1\n");
	fprintf(fp, "cell 2 0 helium 1 -2\n");
	fprintf(fp, "cell 3 0 clad     2 -3\n");
	fprintf(fp, "cell 4 0 water  3 -4\n");
	fprintf(fp, "cell 5 0 outside   4\n");
	fprintf(fp, "\n");
	fprintf(fp, "%% ---------------------------------------------------------------------------- - \n");
	fprintf(fp, "\n");
	fprintf(fp, "%% Material data : \n");
	fprintf(fp, "\n");
	fprintf(fp, "%% ---------------------------------------------------------------------------- - \n");
	fprintf(fp, "\n");
	fprintf(fp, "%% --Fuel :\n");
	fprintf(fp, "\n");
	fprintf(fp, "mat fuel   %E\n", nd_u5 + nd_u8 + nd_n14);
	fprintf(fp, "92235.03c   %E\n", nd_u5);
	fprintf(fp, "92238.03c   %E\n", nd_u8);
	fprintf(fp, "6000.03c   %E\n", nd_n14);
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "%% -- - Carbon buffer layer : \n");
	fprintf(fp, "\n");
	fprintf(fp, "mat buffer   %E     moder grph 6000\n", nd_buf);
	fprintf(fp, "\n");
	fprintf(fp, "6000.03c    %E\n", nd_buf);
	fprintf(fp, "\n");
	fprintf(fp, "%% -- - Pyrolytic carbon layer : \n");
	fprintf(fp, "\n");
	fprintf(fp, "mat PyC     %E      moder grph 6000\n", nd_pyc);
	fprintf(fp, "\n");
	fprintf(fp, "6000.03c    %E\n", nd_pyc);
	fprintf(fp, "\n");
	fprintf(fp, "%% -- - Silicon carbide layer : \n");
	fprintf(fp, "\n");
	fprintf(fp, "mat SiC     %E    moder grph 6000\n", nd_sic * 2);
	fprintf(fp, "\n");
	fprintf(fp, "14028.03c    %.6E\n", 0.9223 * nd_sic);
	fprintf(fp, "14029.03c    %.6E\n", 0.0467 * nd_sic);
	fprintf(fp, "14030.03c    %.6E\n", 0.0310 * nd_sic);
	fprintf(fp, "6000.03c    %E\n", nd_sic);
	fprintf(fp, "\n");
	fprintf(fp, "mat water   %E    moder lwtr 1001\n", nd_h + nd_o);
	fprintf(fp, "8016.03c    %E\n", nd_o);
	fprintf(fp, "1001.03c    %E\n", nd_h);
	fprintf(fp, "\n");
	fprintf(fp, "mat helium   2.6900e-5\n");
	fprintf(fp, "\n");
	fprintf(fp, "2004.03c    1.0\n");
	fprintf(fp, "8016.03c    1E-15\n");
	fprintf(fp, "\n");
	fprintf(fp, "mat clad      0.0970526127000   \n");
	fprintf(fp, "26054.03c    1.896163E-03\n");
	fprintf(fp, "26056.03c    2.976709E-02\n");
	fprintf(fp, "26057.03c    6.876197E-04\n");
	fprintf(fp, "24052.03c    3.235087E-02\n");
	fprintf(fp, "13027.03c    3.235087E-02\n");
	fprintf(fp, "\n");
	fprintf(fp, "%% -- - Thermal scattering data : \n");
	fprintf(fp, "\n");
	fprintf(fp, "therm grph gre7.00t\n");
	fprintf(fp, "therm lwtr lwe7.00t\n");
	fprintf(fp, "%% ---------------------------------------------------------------------------- - \n");
	fprintf(fp, "\n");
	fprintf(fp, "%% Calculation parameters : \n");
	fprintf(fp, "\n");
	fprintf(fp, "%% ---------------------------------------------------------------------------- - \n");
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "%% -- - Geometry plot : \n");
	fprintf(fp, "\n");
	fprintf(fp, "plot 3 500 500\n");
	fprintf(fp, "\n");
	fprintf(fp, "%% -- - Libraries :\n");
	fprintf(fp, "\n");
	fprintf(fp, "set acelib \"/opt/xs/sss_endfb7u.xsdata\"\n");
	fprintf(fp, "set declib \"/opt/xs/sss_endfb7.dec\"\n");
	fprintf(fp, "set nfylib \"/opt/xs/sss_endfb7.nfy\"\n");
	fprintf(fp, "\n");
	fprintf(fp, "%% -- - Boundary conditions : \n");
	fprintf(fp, "\n");
	fprintf(fp, "set bc 2\n");
	fprintf(fp, "\n");
	fprintf(fp, "%% -- - Histories :\n");
	fprintf(fp, "\n");
	fprintf(fp, "set pop 200000 500 50\n");
	fclose(fp);
}

//等效模型输入卡自动生成
void RPT_CARD(double R)
{
	double ND_fuel;
	double ND_U5;
	double ND_U8;
	double ND_N14;
	double ND_Si;
	double ND_C;
	char RPT_filein[500];
	//RPT_filein = (char*)malloc(500);

	ND_U5 = ZL_U5 / (pi * R * R * h);
	ND_U8 = ZL_U8 / (pi * R * R * h);
	ND_N14 = ZL_N14 / (pi * R * R * h);

	ND_Si = (ZL_Si - (r_fcm * r_fcm - R * R) * pi * h * nd_sic) / (pi * R * R * h);
	ND_C = (ZL_C - (r_fcm * r_fcm - R * R) * pi * h * nd_sic) / (pi * R * R * h);
	ND_fuel = ND_U5 + ND_U8 + ND_N14 + ND_Si + ND_C;


	sprintf(RPT_filein, "%f", R);


	fp1 = fopen(RPT_filein, "w");
	fprintf(fp1, "\n");
	fprintf(fp1, "set title \"HTGR Depletion Benchmark, Pebble Bed Model\"\n");
	fprintf(fp1, "\n");
	fprintf(fp1, "%% ---------------------------------------------------------------------------- - \n");
	fprintf(fp1, "\n");
	fprintf(fp1, "%% Geometry:\n");
	fprintf(fp1, "\n");
	fprintf(fp1, "%% ---------------------------------------------------------------------------- - \n");
	fprintf(fp1, "\n");
	fprintf(fp1, "%% -- - Coated particle : \n");
	fprintf(fp1, "\n");
	fprintf(fp1, "pin 1\n");
	fprintf(fp1, "\n");
	fprintf(fp1, "fuel    %f\n", R);
	fprintf(fp1, "SiC     %.4f\n", r_fcm);
	fprintf(fp1, "helium  %.4f\n", r_gap);
	fprintf(fp1, "clad     %.4f\n", r_clad);
	fprintf(fp1, "water\n");
	fprintf(fp1, "\n");
	fprintf(fp1, "%% -- - Universe surrounding particles : \n");
	fprintf(fp1, "\n");
	fprintf(fp1, "%% -- - Pebble\n");
	fprintf(fp1, "surf 4 cuboid %.4f %.4f %.4f %.4f -2 2\n", -1.0 * L, L, -1.0 * L, L);
	fprintf(fp1, "\n");
	fprintf(fp1, "cell 1 0 fill   1 -4\n");
	fprintf(fp1, "cell 5 0 outside   4\n");
	fprintf(fp1, "\n");
	fprintf(fp1, "%% ---------------------------------------------------------------------------- - \n");
	fprintf(fp1, "\n");
	fprintf(fp1, "%% Material data : \n");
	fprintf(fp1, "\n");
	fprintf(fp1, "%% ---------------------------------------------------------------------------- - \n");
	fprintf(fp1, "\n");
	fprintf(fp1, "%% --Fuel :\n");
	fprintf(fp1, "\n");

	fprintf(fp1, "mat fuel    %E  moder grph 6000\n", ND_fuel);
	fprintf(fp1, "92235.03c   %E\n", ND_U5);
	fprintf(fp1, "92238.03c   %E\n", ND_U8);
	fprintf(fp1, "14028.03c    %.6E\n", 0.9223 * nd_sic);
	fprintf(fp1, "14029.03c    %.6E\n", 0.0467 * nd_sic);
	fprintf(fp1, "14030.03c    %.6E\n", 0.0310 * nd_sic);
	fprintf(fp1, " 6000.03c   %E\n", ND_C+ND_N14);
	fprintf(fp1, "\n");
	fprintf(fp1, "mat SiC     %E    moder grph 6000\n", nd_sic * 2);
	fprintf(fp1, "\n");
	fprintf(fp1, "14028.03c    %.6E\n", 0.9223 * ND_Si);
	fprintf(fp1, "14029.03c    %.6E\n", 0.0467 * ND_Si);
	fprintf(fp1, "14030.03c    %.6E\n", 0.0310 * ND_Si);
	fprintf(fp1, "6000.03c    %E\n", nd_sic);
	fprintf(fp1, "\n");
	fprintf(fp1, "\n");
	fprintf(fp1, "mat water   %E    moder lwtr 1001\n", nd_h + nd_o);
	fprintf(fp1, "8016.03c    %E\n", nd_o);
	fprintf(fp1, "1001.03c    %E\n", nd_h);
	fprintf(fp1, "\n");
	fprintf(fp1, "mat helium   2.6900e-5\n");
	fprintf(fp1, "\n");
	fprintf(fp1, "2004.03c    1.0\n");
	fprintf(fp1, "8016.03c    1E-15\n");
	fprintf(fp1, "\n");
	fprintf(fp1, "mat clad      0.0970526127000   \n");
	fprintf(fp1, "\n");
	fprintf(fp1, "26054.03c    1.896163E-03\n");
	fprintf(fp1, "26056.03c    2.976709E-02\n");
	fprintf(fp1, "26057.03c    6.876197E-04\n");
	fprintf(fp1, "24052.03c    3.235087E-02\n");
	fprintf(fp1, "13027.03c    3.235087E-02\n");
	fprintf(fp1, "\n");
	fprintf(fp1, "%% -- - Thermal scattering data : \n");
	fprintf(fp1, "\n");
	fprintf(fp1, " therm grph gre7.00t\n");
	fprintf(fp1, "therm lwtr lwe7.00t\n");
	fprintf(fp1, "%% ---------------------------------------------------------------------------- - \n");
	fprintf(fp1, "\n");
	fprintf(fp1, "%% Calculation parameters : \n");
	fprintf(fp1, "\n");
	fprintf(fp1, "%% ---------------------------------------------------------------------------- - \n");
	fprintf(fp1, "\n");
	fprintf(fp1, "\n");
	fprintf(fp1, "%% -- - Geometry plot : \n");
	fprintf(fp1, "\n");
	fprintf(fp1, "plot 3 500 500\n");
	fprintf(fp1, "\n");
	fprintf(fp1, "%% -- - Libraries :\n");
	fprintf(fp1, "\n");
	fprintf(fp1, "set acelib \"/opt/xs/sss_endfb7u.xsdata\"\n");
	fprintf(fp1, "set declib \"/opt/xs/sss_endfb7.dec\"\n");
	fprintf(fp1, "set nfylib \"/opt/xs/sss_endfb7.nfy\"\n");
	fprintf(fp1, "\n");
	fprintf(fp1, "%% -- - Energy grid : \n");
	fprintf(fp1, "\n");
	fprintf(fp1, "\n");
	fprintf(fp1, "%% -- - Boundary conditions : \n");
	fprintf(fp1, "\n");
	fprintf(fp1, "set bc 2\n");
	fprintf(fp1, "\n");
	fprintf(fp1, "%% -- - Histories :\n");
	fprintf(fp1, "\n");
	fprintf(fp1, "set pop 200000 500 50\n");

	fclose(fp1);

}


//读取k
void ReadFCMK()
{

	int i;
	char none[500];
	char symbol[] = "IMP_KEFF";
	fp2 = fopen("FCM_res.m", "r");
	for (i = 0; i < inf; i++)
	{
		fscanf(fp2, "%s", none);
		if (strcmp(none, symbol) == 0)
		{
			//IMP_KEFF                  (idx, [1:   2]) = [  7.32266E-01 0.00022 ];

			fscanf(fp2, "%s", none);//IMP_KEFF
			fscanf(fp2, "%s", none);//(idx,
			fscanf(fp2, "%s", none);// [1:
			fscanf(fp2, "%s", none);//2])
			fscanf(fp2, "%s", none);//=
			fscanf(fp2, "%lf", &K_fcm);//[ 
			fscanf(fp2, "%lf", &bias_fcm);
		}
	}


}

//读取k
void ReadRPTK(double& K_rpt, double& bias_rpt, double R)
{

	int i;
	char none[500];
	char symbol[] = "IMP_KEFF";
	char RPT_filein[500];
	char RPT_fileout[500];
	sprintf(RPT_filein, "%f", R);
	sprintf(RPT_fileout, "%s_res.m", RPT_filein);
	printf("%s\n", RPT_fileout);
	fp2 = fopen(RPT_fileout, "r");
	for (i = 0; i < inf; i++)
	{
		fscanf(fp2, "%s", none);
		if (strcmp(none, symbol) == 0)
		{
			fscanf(fp2, "%s", none);
			fscanf(fp2, "%s", none);
			fscanf(fp2, "%s", none);
			fscanf(fp2, "%s", none);
			fscanf(fp2, "%s", none);
			fscanf(fp2, "%lf", &K_rpt);
			printf("k=%f\n", K_rpt);
			fscanf(fp2, "%lf", &bias_rpt);
		}
	}


}

int main()
{
	FILE* tp;
	tp = fopen("record.m", "w");
	fprintf(tp, "Record the RPT radius and kinf\n\n");
	char run_command[100];
	char RPT_filein[500];
	int i, j, k;
	double Rmin, Rmax;
	Read_para();
	Rmin = 0.6324555 * r_fcm;      //半径下限（这么说就是40%的填充率咯）
	Rmax = r_fcm;                //半径上限

	CalZL();
	FCM_CARD();
	printf("FCM CARD OK\n");
	system("sss2 -omp 38 FCM");  //运行serpent
	ReadFCMK();
	printf("K_fcm=%f\n", K_fcm);


	R2 = Rmax;
	RPT_CARD(R2);
	sprintf(RPT_filein, "%f", R2);
	sprintf(run_command, "sss2 -omp 38 %s", RPT_filein);
	system(run_command);//运行serpent
	ReadRPTK(K_rpt2, bias_rpt2, R2);
	printf("K2=%f  Kfcm=%f\n", K_rpt2, K_fcm);
	deltK2 = (K_rpt2 - K_fcm) * 100000;

	R1 = Rmin;
	RPT_CARD(R1);
	sprintf(RPT_filein, "%f", R1);
	sprintf(run_command, "sss2 -omp 38 %s", RPT_filein);
	system(run_command); //运行serpent
	ReadRPTK(K_rpt1, bias_rpt1, R1);
	deltK1 = (K_rpt1 - K_fcm) * 100000;



	//搜索半径
	for (i = 0; i < cyc_max; i++)
	{

		printf("K1=%f  K2=%f  Kfcm=%f\n", K_rpt1, K_rpt2, K_fcm);
		R3 = (R2 - R1) * (K_fcm - K_rpt1) / (K_rpt2 - K_rpt1) + R1;

		printf("R3=%f\n", R3);
		RPT_CARD(R3);
		sprintf(RPT_filein, "%f", R3);
		//sprintf(run_command, "sss2 -omp 40 %s", RPT_filein);
		sprintf(run_command, "sss2 -omp 38 %s", RPT_filein);
		system(run_command);  //运行serpent
		ReadRPTK(K_rpt3, bias_rpt3, R3);
		deltK3 = (K_rpt3 - K_fcm) * 100000;

		printf("deltK3=%.2f\n", deltK3);
		fprintf(tp, "i=%2d  R_rpt=%f  k_fcm=%f  k_rpt=%f  bias_fcm=%f  bias_rpt=%f\n", i, R3, K_fcm, K_rpt3, bias_fcm, bias_rpt3);
		if ((deltK3 >= (-1.0 * pcm)) && (deltK3 <= pcm))
		{
			i = 100;
		}

		if (deltK3 > pcm)
		{
			R1 = R3;
			K_rpt1 = K_rpt3;
		}
		else if (deltK3 < (-1.0 * pcm))
		{
			R2 = R3;
			K_rpt2 = K_rpt3;
		}


	}
	fclose(tp);
	printf("Finish searching RPT Radius!\n");
	printf("The Radius of RPT method : %f\n", R3);
	printf("The kinf of RPT method : %f  %f\n", K_rpt3, bias_rpt3);
	printf("The kinf of FCM fuel :   %f  %f\n", K_fcm, bias_fcm);

}