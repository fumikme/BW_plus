/** �X�V���� **/
/********* 2019/04/09  **********///											 �������������؂�܂�
/* FEM_LPmode_VCSELinput���A�R�A���a10um�̃X�e�b�v�^���z�t�@�C�o�̕��z�Ƃ��Čv�Z����FEM_result����Z�o���Ă���_�ɒ��� */

/*******************************************************************************************************/
/*******************************                              ******************************************/
/*******************************         �֐S���̕���         ******************************************/
/*******************************                              ******************************************/
/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
//                                                                                                     //
/*-*-                     ���W���[�������Ă��ꂼ��10�s���炢�ɂ��Ă����ƁC                        -*-*-*/
/*-*-             �������鎞�ɂ����������������P�ʂ̃��W�b�N�Ɏv�l���W���ł���̂ŁA              -*-*-*/
/*-*-                              �R�[�h�������₷���Ȃ�܂��B                                   -*-*-*/
//                                                                                                     //
/*******************************************************************************************************/

#define _USE_MATH_DEFINES
#define NR_END 1
#define FREE_ARG char*
#pragma warning(disable:4996) //fopen�֐�������
#pragma warning(disable:0266)

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;


/* ���萔�̒�` */
#define PI 3.141592653 // �~����
#define C 2.99792458e8 // ����
#define Epsilon0 8.8541878e-12 // �^�󒆂̗U�d��
#define Mu0 PI*4.0e-7 // �^�󒆂̓�����

#define num_mode_vin 6

/* �֐����v���O�����̓ǂݍ��� *//* �֐����v���O�����̓ǂݍ��� */
// ���ܗ��g������
double dndl(double lamda, double n_lamda, int mater);
// �W���s��v�f�v�Z�֐�
void S_matrix(double *a, double *b, double *q, int m, int n, double v, double w, double D);
// �����R���X�L�[��������щ����R���X�L�[�����ɂ�鋁��
void mcholesky(double *a, double *b, double *ML, double *MD, int m, int n);
void mcholesky_sol(double *ML, double *MD, double *R, int m, int n);
// �t�ׂ���@�̏����x�N�g���v�Z�E���x�N�g���K�i��
void R0(double *MD, double *R, int m, int n);
void R_norm(double *R, int n);
// �ŗL�l�v�Z
double Eigen(double *R, double *a, double *b, int n, int m);
//  �Q�x���v�Z�p�֐�
double dbdk_bunbo(double *R, double D, double w, int m, int n);
double dbdk_bunshi(double *R, double *qg, double D, double w, int m, int n);
// ��ݍ��ݐϕ�
double *convolution(int xmax, double *P, int LT, double *pulse);

/* ���w�֐��E�������z��̓ǂݍ��� */
// ��1��ό`Bessel�֐� In(x)�E��2��ό`Bessel�֐� Kn(x)
double bessi0(double x), bessi1(double x);
double bessk0(double x), bessk1(double x), bessk(int n, double x);
// 1�����z��i�����^�j
int  *dintvector(int i, int j);
void free_dintvector(int *a, int i);
void init_intvector(int *a, int nr1, int nr2);
// 1�����z��i�����^�j
double *drealvector(int i, int j);
void free_drealvector(double *a, int i);
void init_realvector(double *a, int nr1, int nr2);
// 2�����z��ifloat�^�j
float **matrix(long nrl, long nrh, long ncl, long nch);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void init_matrix(float **a, int nr1, int nr2, int nl1, int nl2);
// 2�����z��idouble�^�j
double  **dmatrix(int nr1, int nr2, int nl1, int nl2);
void free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2);
void init_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2);
// 3�����z��
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void init_d3tensor(double ***a, int nr1, int nr2, int nl1, int nl2, int np1, int np2);

/* �G���[�o�͊֐��̓ǂݍ��� */
void nrerror(char error_text[]);

/* 0. �p�����[�^�̒�` */
// m: ���[�h����( TE&TM ~ 1, EH ~ n+1, HE ~ n-1), l: ���a�������[�h����, n: ���ʊp���[�h����
// mater: �ޗ�ID, profile: ���ܗ����z���͕��@ID (0: �ׂ��敪�z, 1: ���蕪�z����)
// NLP: LP���[�h��, NLP0: 0�����ʊp���[�hLP���[�h��, Ntotal: �S���[�h��
// N: ���a���Wr�̃R�A����ԃX�e�b�v��, Nclad: ��͗̈��ԃX�e�b�v��, Nbeta: �`���萔�̕�����
// Nwkb: WKB�ߎ��ɂ�苁�߂�LP���[�h��, Ntmin: �ŏ��Q�x�����ԃX�e�b�v��, Ntmax: �ő�Q�x�����ԃX�e�b�v��
// lamda: �g��, k: �g��, omega: �p���g��, delta: ����ܗ�, NA: �J����, aa: ���a���W�K�i���T�C�Y
// A: �R�A���a, AA: ��͗̈攼�a, g: �ׂ��敪�z�w��, n1: �R�A���S�̋��ܗ�, n0: �N���b�h�̋��ܗ�
// L: �t�@�C�o��, Tv: �C���p���X�������Ԏ����ݕ�, dr: ���a���W���ݕ�,
// v: �K�i�����g��, w: �K�i���`���萔, D: �K�i���R�A�a
// tau: �Q���x, beta: �`���萔��, dbeta: �`���萔���ݕ�
// a: �W���s��S�̕��Ίp�v�f�i�[�z��, b: �W���s��S�̑Ίp�v�f�i�[�z��
// GI: �R�A�����ܗ�, q: �K�i���R�A�����ܗ�, qg: �R�A�����ܗ����U�p�����[�^
// R: �K�i���������d�ꐬ��, Rb: �t�ׂ���@�p����q�z��, Rinf: ��͗̈�
// ML: LDU�����W���s���L�s�񕛑Ίp�v�f, MD: LDU�����W���s���D�s��Ίp�v�f
// de, df: �t�ׂ���@�ɂ�����������ʃp�����[�^, eig: �ŗL�l
// taumin: �ŏ��P�ʌQ�x��, taumax: �ő�P�ʌQ�x��
// fmin: �]���ŏ����g��, fmax: �]���ő���g��, dfrq: �]�����g������\
// Hw: �`�B�֐��v�f����q�ϐ� (ReHw: ����, ImHw: ����), bw: -3dB�ш敝
// GI: ���ܗ����z, q: ���ܗ��g������, qg: n*d(kn)/dk
// Rlp: �SLP���[�h�̃��[�h�t�B�[���h, Mbeta: �e���[�h�̓`���萔
// r0: �I�t�Z�b�g�i��U�ʒu�j�Cw0: �K�E�V�A���r�[�����a�Cdx: �d�Ȃ�ϕ��̋�ԕ���\, Nxy: ��͗̈�̕�����


int main(void)
{
	FILE   *fp, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7, *fp8, *fr, *fpulse, *fopulse, *fmpd;
	char   fppath[128], fp2path[128], fp3path[128], fp4path[128], fp5path[128], fp6path[128], fp7path[128];
	char   fp8path[128], frpath[128], fpulsepath[128], fopulsepath[128], fmpdpath[128];
	          fppath[0]   =	 fp2path[0] = fp3path[0] = fp4path[0] = fp5path[0] = fp6path[0] = fp7path[0] =
	          fp8path[0]  =  frpath[0] = fpulsepath[0] = fopulsepath[0] = fmpdpath[0] = '\0' ;
	int    myu, x, y, i, j, jmax, count;
	int    m, l, NLP, NLP0, Ntotal, N, Nclad, Nbeta, mater, profile, Nwkb, Ntmin, Ntmax, Nf, launch;
	int    Nxy, xmax, LT, xL;
	double lamda, k, omega, A, AA, g, n0, n1, dr, L, Tv;
	double r0, w0, dx, dy;
	double delta, NA, aa, v, w, D;
	double tau, beta, dbeta, bb, eps1, eps2, sum, Rinf, de, df, eig;
	double taumin, taumax, fmin, fmax, dfrq, Hw, ReHw, ImHw, bw;					//	int    nr;
	double yy, Rxy, Rinxy;
	double Em_evin, Em_odin, Em_ev, Em_od, Am_evev, Am_evod, Am_odev, Am_odod;
	double ncav, noxi, Avin;
	double *GI, *pulse, *opulse, *q, *qg, *R, *R2, *Rb, *a, *b, *ML, *MD, *Mtau, *P, *M, *Mbeta, *MPD;
	double *Mbetain;
	double **Rlp, **MPD2d, **Rinlp;;
	double tauin, betain, Rinfin, eigin, betain_devided_by_k, data;
	double win, xx1, xx2, rr1, rr2;
	int    *modem, *modemin, *model, *modelin;
	int    Nvin, couple, min, lin, nrr1, nrr2, max;
	int    OFFres, OFFrange;

	char   trash[65536] = "\0";
	char   directory[128] = "FILE_for_IO";

	/** 1. ���̓t�@�C���̓ǂݍ��� **/
	if ((fp = fopen("FEM_input_kob.csv", "r")) != NULL) {
		char s1[256], s2[256];
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &A);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &AA);
		fscanf(fp, "%[^,], %[^,], %d\n", s1, s2, &profile);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &g);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &n1);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &n0);
		fscanf(fp, "%[^,], %[^,], %d\n", s1, s2, &launch);
		fscanf(fp, "%[^,], %[^,], %d\n", s1, s2, &couple);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &r0);
		fscanf(fp, "%[^,], %[^,], %d\n", s1, s2, &OFFrange);
		fscanf(fp, "%[^,], %[^,], %d\n", s1, s2, &OFFres);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &Avin);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &w0);			// �K�E�X���z�i�V���O�����[�h���[�U�j��p���ė�U
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &dr);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &dx);
		fscanf(fp, "%[^,], %[^,], %d\n", s1, s2, &mater);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &lamda);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &L);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &Tv);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &fmin);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &fmax);
		fscanf(fp, "%[^,], %[^,], %d\n", s1, s2, &Nf);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &aa);
		fscanf(fp, "%[^,], %[^,], %d\n", s1, s2, &Nbeta);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &eps1);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &eps2);
		fscanf(fp, "%[^,], %[^,], %d\n", s1, s2, &jmax);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &noxi);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &ncav);
		// fscanf(fp, "%[^,]\n", s1);
	}
	else { printf(" U cannot open the file !\n"); exit(EXIT_FAILURE); }

	N = (int)(1000.0*A / dr); Nclad = (int)(1000.0*(AA - A) / dr);				//-----------Nclad�ɒ���
	Nvin = (int)(1000.0*Avin / dr);


	fmin = fmin * 1.0e-3; fmax = fmax * 1.0e-3; dfrq = (fmax - fmin) / (double)Nf;

	/************************   ���������ɏd�v   ****************************/
	Nxy = (int)(3 * 1000 * A / dx); //���胂�[�h��U���̉�͗̈�ɑΉ����镪���� Nxy*dx=��͗̈�			
	/**************************************************************************/

	lamda = lamda * 1.0e-9;  A = A * 1.0e-6; dr = dr * 1.0e-9; AA = AA * 1.0e-6;
	Avin = Avin * 1.0e-6;
	r0 = r0 * 1.0e-6; w0 = w0 * 1.0e-6, dx = dx * 1.0e-9; dy = dx;	// �P�ʕϊ��iTHz(1/ps)�j			
	fclose(fp);

	/* ����v���t�@�C���ǂݍ��� */
	GI = drealvector(0, N); init_realvector(GI, 0, N);
	sprintf(frpath, "%s/profile.csv", directory);
	if (profile == 1) {
		if ((fr = fopen(frpath, "r")) != NULL) {
			for (j = 0; j <= N; j++) { fscanf(fr, "%lf,", &GI[j]); }
		}
		else { printf(" U cannot open the file !\n"); exit(EXIT_FAILURE); }
		n1 = GI[0]; n0 = GI[N]; fclose(fr);
	}

	/* ���˔g�`�ǂݍ��� */
	pulse = drealvector(0, 1999); init_realvector(pulse, 0, 1999); //���˔g�`��2000�X�e�b�v�œ���
	sprintf(fpulsepath, "%s/input pulse.csv", directory);
	if ((fpulse = fopen(fpulsepath, "r")) != NULL) {
		for (j = 0; j <= 1999; j++) { fscanf(fpulse, "%lf,", &pulse[j]); }
	}
	else { printf(" U cannot open the file !\n"); exit(EXIT_FAILURE); }
	fclose(fpulse);


	/** 2. �e��萔�̐ݒ� **/
	k = 2.0*PI / lamda; omega = 2.0*PI*C / lamda;
	delta = (n1*n1 - n0 * n0) / (2.0*n1*n1); NA = sqrt(n1*n1 - n0 * n0);
	v = k * aa*n1*sqrt(2.0*delta); D = A / aa;

	Nwkb = int((1.0 / 4.0)*(g / (g + 2.0))*(k*k)*(n1*n1)*delta*(A*A));
	bb = -1; dbeta = k * (n1 - n0) / (double)Nbeta;//������
	NLP = NLP0 = Ntotal = 0;//������

							/* �z��̋L���̈�m�ۂ���я����� */
	q = drealvector(0, N); init_realvector(q, 0, N);
	qg = drealvector(0, N); init_realvector(qg, 0, N);
	R = drealvector(0, N); init_realvector(R, 0, N);
	R2 = drealvector(0, N); init_realvector(R2, 0, N);
	Rb = drealvector(0, N); init_realvector(Rb, 0, N);
	a = drealvector(0, N); init_realvector(a, 0, N);
	b = drealvector(0, N); init_realvector(b, 0, N);
	ML = drealvector(0, N); init_realvector(ML, 0, N);
	MD = drealvector(0, N); init_realvector(MD, 0, N);
	M = drealvector(0, Nf); init_realvector(M, 0, Nf);
	Mtau = drealvector(0, 10 * Nwkb); init_realvector(Mtau, 0, 10 * Nwkb);
	modem = dintvector(0, 10 * Nwkb); init_intvector(modem, 0, 10 * Nwkb);
	model = dintvector(0, 10 * Nwkb); init_intvector(model, 0, 10 * Nwkb);
	modelin = dintvector(0, 10 * Nwkb); init_intvector(modelin, 0, 10 * Nwkb);
	Mbeta = drealvector(0, 10 * Nwkb); init_realvector(Mbeta, 0, 10 * Nwkb);
	MPD = drealvector(0, 10 * Nwkb); init_realvector(MPD, 0, 10 * Nwkb);

	#define sizeofRLP 500
	Rlp = dmatrix(0, sizeofRLP, 0, sizeofRLP);
	for (int i = 0; i <= sizeofRLP; i++) {
		for (int j = 0; j <= sizeofRLP; j++) { Rlp[i][j] = 0.0; }
	}

	/* ���ܗ����z�̐ݒ� */
	/* �ׂ���v���t�@�C��*/
	if (profile == 0) {
		for (j = 0; j <= N; j++) { GI[j] = n1 * sqrt(1.0 - 2.0*delta*pow(((double)j / (double)N), g)); }
	}
	for (j = 0; j <= N; j++) { q[j] = (GI[j] * GI[j] - n0 * n0) / (n1*n1 - n0 * n0); }
	for (j = 0; j <= N; j++) { qg[j] = GI[j] * GI[j] - (lamda*GI[j] * dndl(lamda*1.0e9, GI[j], mater)) / (1 - (lamda / GI[j])*dndl(lamda*1.0e9, GI[j], mater)); }


	/** 3. �]�������̏o�� **/
	sprintf(fp2path, "%s/FEM_setting.csv", directory);
	if ((fp2 = fopen(fp2path, "w")) != NULL) {
		fprintf(fp2, "Material, mater, %d\n", mater);
		fprintf(fp2, "Wavelength, ��, %lf,nm\n", lamda*1e9);
		fprintf(fp2, "Fiber length, L, %lf,m\n", L);
		fprintf(fp2, "Core radius, A, %lf,��m\n", A*1e6);
		fprintf(fp2, "Analysis region in radial axis, AA,%lf,��m\n", AA*1e6);
		fprintf(fp2, "Index exponent, g, %lf\n", g);
		fprintf(fp2, "Refractive index at the core center, n1,%lf\n", n1);
		fprintf(fp2, "Refractive index in the cladding, n0,%lf\n", n0);
		fprintf(fp2, "Relative refractive index, ��,%lf\n", delta);
		fprintf(fp2, "Numerical aperture, NA, %lf\n", NA);
		fprintf(fp2, "Step size of the elements, dr, %lf, nm\n", dr*1e9);
		fprintf(fp2, "Step size of the pulse waveform, Tv, %lf, ps\n", Tv);
		fprintf(fp2, "Minimum evaluated frequency,fminx,%e,GHz\n", fmin*1.0e3);
		fprintf(fp2, "Maximum evaluated frequency,fmax,%e,GHz\n", fmax*1.0e3);
		fprintf(fp2, "Step size of evaluated frequency,dfrq,%e,GHz\n", dfrq*1.0e3);
		fprintf(fp2, "Step size of propagation constant, d��,%lf\n", dbeta);
		fprintf(fp2, "Partition number of propagation constant, N��,%d\n", Nbeta);
		fprintf(fp2, "Partition number of fiber core radius, N,%d\n", N);
		fprintf(fp2, "Partition number of fiber cladding, Nclad,%d\n", Nclad);
		fprintf(fp2, "Maximum allowable error for convergence solution vector, eps1,%lf\n", eps1);
		fprintf(fp2, "Maximum allowable error of zero eigen value, eps2,%lf\n", eps2);
		fprintf(fp2, "Maximum number of iterations in inverse power method, jmax,%d\n", jmax);
		printf("Material no.: %d\n", mater);
		printf("Wavelength ��: %lfnm\n", lamda*1.0e9);
		printf("Fiber length L: %lfm\n", L);
		printf("Core radius A: %lf��m\n", A*1.0e6);
		printf("Analysis region AA: %lf ��m\n", AA*1e6);
		printf("Index exponent g: %lf\n", g);
		printf("Refractive index at the core center n1: %lf\n", n1);
		printf("Refractive index in the cladding n0: %lf\n", n0);
		printf("Relative refractive index ��: %lf\n", delta);
		printf("Numerical aperture NA: %lf\n", NA);
		printf("Step size of the elements dr: %lfnm\n", dr*1.0e9);
		printf("Partition number of fiber core radius N: %d\n", N);
		printf("Partition number of fiber cladding Nclad: %d\n", Nclad);
		printf("Step size of propagation constants d��: %lf\n", dbeta);
		printf("Partition number of propagation constants N��: %d\n", Nbeta);
		printf("Maximum allowable error for convergence solution vector eps1: %lf\n", eps1);
		printf("Maximum allowable error of zero eigen value eps2: %lf\n", eps2);
		printf("Maximum number of iterations in inverse power method, jmax,%d\n", jmax);
	}
	else { printf(" U cannot open the file !\n"); exit(EXIT_FAILURE); }

	sprintf(fp3path, "%s/FEM_profile.csv", directory);
	if ((fp3 = fopen(fp3path, "w")) != NULL) {
		fprintf(fp3, "r [��m], n, q, qg\n");
		for (j = 0; j <= N; j++) {
			fprintf(fp3, "%lf, %lf, %lf, %lf\n", A*1.0e6* (double)j / (double)N, GI[j], q[j], qg[j]);
		}
	}
	else { printf(" U cannot open the file !\n"); exit(EXIT_FAILURE); }
	fclose(fp3);

	if ((fp4 = fopen("FEM_result.csv", "w")) != NULL) {
		printf("m, l, tau, beta, ne, eig, R[N]\n");
		fprintf(fp4, "m,l,tau,beta,ne,R_infinite,eig,r=0_um,");
		for (j = 1; j <= N; j++) { fprintf(fp4, "%lf,", (double)j*dr*1.0e6); } fprintf(fp4, "\n");
	}
	else { printf(" U cannot open the file !\n"); exit(EXIT_FAILURE); }

	sprintf(fp5path, "%s/FEM_impulse_response.csv", directory);
	if ((fp5 = fopen(fp5path, "w")) != NULL) { ; }
	else { printf(" U cannot open the file !\n"); exit(EXIT_FAILURE); }

	sprintf(fp6path, "%s/FEM_frequency-respnse.csv", directory);
	if ((fp6 = fopen(fp6path, "w")) != NULL) {
		for (y = 0; y <= Nf; y++) { fprintf(fp6, "%lf,", (fmin*1.0e3) + (double)y*(dfrq*1.0e3)); } fprintf(fp6, "\n");
	}
	else { printf("The file cannot be opened !\n"); exit(EXIT_FAILURE); }

	sprintf(fopulsepath, "%s/FEM_output pulse.csv", directory);
	if ((fopulse = fopen(fopulsepath, "w")) != NULL) { ; }
	else { printf(" U cannot open the file !\n"); exit(EXIT_FAILURE); }

	sprintf(fmpdpath, "%s/FEM_mode power distribution.csv", directory);
	if ((fmpd = fopen(fmpdpath, "w")) != NULL) { ; }
	else { printf(" U cannot open the file !\n"); exit(EXIT_FAILURE); }
	

	/** 4. �ŗL�x�N�g������ьŗL�l�̌v�Z **/
	for (m = 0; ; m++) {
		l = 1; beta = k * n1;

		/****************************************************************************************/
		for (; ; ) {
			/* �W���s��v�Z����щ����R���X�L�[���� */
			w = aa * sqrt((beta*beta) - (k*k)*(n0*n0));
			S_matrix(a, b, q, m, N, v, w, D);
			mcholesky(a, b, ML, MD, m, N);
			/* �����x�N�g�� R0 �̕t�^ */
			R0(MD, R, m, N);
			/*  �A���ꎟ������ SR=(LDL)R=bR �̔����]�� */
			for (j = 0; ; j++) {
				for (i = 0; i <= N; i++) { Rb[i] = R[i]; }
				mcholesky_sol(ML, MD, R, m, N);
				R_norm(R, N);
				/* �������� */
				de = 0, df = 0;
				for (i = 0; i <= N; i++) {
					de = de + (Rb[i] - R[i])*(Rb[i] - R[i]);
					df = df + (Rb[i] + R[i])*(Rb[i] + R[i]);
				}

				if (de < eps1 || df < eps1) break;
				if (j >= jmax) goto nextin;
				// �@ R��Rb�̐�����de��0�ɑQ�߂���Ύ����ibreak�j�D
				// �A R��Rb�̐����adf��0�ɑQ�߂���Β��~�ibreak�j�D
				// �B �����񐔂�����l jmax �𒴂��������ύX���čČv�Z�igo to next�j�D
			}

			/* �ŗL�l�̌v�Z */
			eig = Eigen(R, a, b, N, m);

			/* �ŗL�l�̑Ó����]�� */
			// �u�����ŗL�l eig ���O��l bb �Ɠ��l�v�ł���΃���ς��ď��߂���Čv�Z
			if (eig == bb) {
				dbeta = k * (n1 - n0) / (double)Nbeta;
				beta = beta - 1.0*dbeta;
				continue;
			}

			//  �@�u0< �����ŗL�l eig < 0.00001�v�ł���Η�ŗL�l�Ƃ��č̗p
			if (0.0 < eig && eig < eps2) {
				/* �������d�ꐬ��R�̋K�i�� �i�p���[��1W�Ƃ���j*/
				sum = 0.0;
				for (j = 0; j < N; j++) { sum = sum + R[j] * R[j] * (j*dr)*dr; }
				for (j = 0; j < Nclad; j++) { sum = sum + R[N] * (bessk(m, w*(j*dr + A)) / bessk(m, w*A))*R[N] * (bessk(m, w*(j*dr + A)) / bessk(m, w*A))*(j*dr + A)*dr; }
				for (j = 0; j <= N; j++) { R2[j] = R[j] * sqrt((omega*Mu0) / (PI*beta*sum)); }
				tau = Mtau[NLP] = (1.0 / (C*1.0e-12))*(k / beta) * dbdk_bunshi(R, qg, D, w, m, N) / dbdk_bunbo(R, D, w, m, N); // = (1/c)*(d��/dk) [ps/m]
				modem[NLP] = m;			model[NLP] = l;
				Mbeta[NLP] = beta;

				/* �v�Z���ʂ̏o�� */
				Rinf = R[N] * (bessk(m, w*(Nclad*dr + A)) / bessk(m, w*A));
				fprintf(fp4, "%d, %d, %lf, %lf, %lf, %lf, %lf,", m, l, tau, beta, beta / k, Rinf, eig);
				for (j = 0; j <= N; j++) {
					fprintf(fp4, "%lf,", R2[j]);
					Rlp[NLP][j] = R2[j];
				}

				fprintf(fp4, "\n");
				printf("%d, %d, %lf, %lf, %lf, %lf, %lf\n", m, l, tau, beta, beta / k, eig, R[N]);
				dbeta = k * (n1 - n0) / (double)Nbeta;
				bb = eig;
				l = l + 1;
				NLP = NLP + 1; if (m == 0) { NLP0 = NLP0 + 1; }
				count = 0;
			}

			//  �A�u0 < �����ŗL�l eig�v���u-1 < �O��l bb < 0�v�ł����
			if (eig > 0.0 && bb < 0.0 && (fabs(bb) < 1.0)) {
				beta = beta + dbeta;
				dbeta = dbeta / 2.0;
				count = count + 1;
			}

			// �B ���̑�
			else { bb = eig; count = 0; }
			if (count > 3000) {
				dbeta = k * (n1 - n0) / (double)Nbeta;
				beta = beta - dbeta;
			}

		nextin:
			beta = beta - dbeta;
			if (beta < k*n0)	break;
		}
		/****************************************************************************************/
		if (l == 1) { Ntotal = (NLP0) * 2 + (NLP - NLP0) * 4; break; }
	}		// break�ƂȂ�̂́Cm���ō������ɓ��B�����Ƃ�
	printf("�ŗL�l�v�Z����\n");

	/* ��U�����ݒ�iMPD�z��̎Z�o�j*/
	// OFL condition
	if (launch == 0) {
		for (m = 0; m < NLP; m++) {
			if (modem[m] == 0) { MPD[m] = 2.0; }
			else { MPD[m] = 4.0; }
		}
	}	  //�k�ސ��̍l��


		  // RML condition																		
#if 1

	if (launch == 1) {
		double Amev, Amod, Emev, Emod, Ein, xx, rr; //��U�����̂Ƃ��ŏo�Ă���
		int nr;
		for (m = 0; m < NLP; m++) {
			Amev = Amod = 0.0;
			for (i = 0; i < Nxy; i++) {
				xx = (r0 - ((double)Nxy*dx / 2.0)) + ((double)i*dx);
				for (j = 0; j < Nxy; j++) {
					yy = (-((double)Nxy*dx / 2.0)) + ((double)j*dy);
					rr = sqrt(xx*xx + yy * yy);
					nr = (int)(rr / dr);
					if (nr == 0) { Rxy = Rlp[m][0] + (Rlp[m][1] - Rlp[m][0]) * (rr / dr); }	
					else if (nr <= N) { Rxy = Rlp[m][nr] + (Rlp[m][nr + 1] - Rlp[m][nr]) * ((rr - (double)nr*dr) / dr); } //�R�A
					else if (nr > N) {
						w = aa * sqrt(Mbeta[m] * Mbeta[m] - k * k*n0*n0);
						Rxy = Rlp[m][N] * (bessk(m, w*rr) / bessk(m, w*A));
					}			//�N���b�h
					Emev = Rxy * cos((double)(modem[m]) * atan2(yy, xx));					
					Emod = Rxy * sin((double)(modem[m]) * atan2(yy, xx));					
					Ein = exp(-((xx - r0)*(xx - r0) + yy * yy) / (w0*w0));					
					Amev = Amev + Emev * Ein*dx*dy;											
					Amod = Amod + Emod * Ein*dx*dy;
				}
			}

			printf("%d ", m);
			Amev = Amev * Amev / (((2.0*omega*Mu0) / Mbeta[m]) * (PI*w0*w0 / 2.0));
			Amod = Amod * Amod / (((2.0*omega*Mu0) / Mbeta[m]) * (PI*w0*w0 / 2.0));			
			MPD[m] = 100.0*(Amev + Amod);
		}
		for (i = 0; i <= NLP - 1; i++) { fprintf(fmpd, "%lf\n", MPD[i]); }
	}
#endif

	printf("%d\n", NLP);


	// RML condition for VCSEL ***********************************************************************
	min = lin = max = nrr1 = nrr2 = 0;
	tauin = betain = Rinfin = eigin = betain_devided_by_k = data = xx1 = xx2 = rr1 = rr2 = 0.0;
	Rinlp = dmatrix(0, (int)(2 * Avin / dr), 0, (int)(2 * Avin / dr));
	MPD2d = dmatrix(0, NLP, 0, NLP);



	modemin = dintvector(0, 10 * Nwkb); init_intvector(modemin, 0, 10 * Nwkb);
	Mbetain = drealvector(0, 10 * Nwkb); init_realvector(Mbetain, 0, 10 * Nwkb);
	
	for (i = 0; i <= (int)(2 * Avin / dr); i++) { for (j = 0; j <= (int)(2 * Avin / dr); j++) { 
		Rinlp[i][j] = 0.0; 		}}
	for (i = 0; i <= NLP; i++) { for (j = 0; j <= NLP; j++) { MPD2d[i][j] = 0.0; } }
	double cef_evev, cef_evod, cef_odev, cef_odod;

	/*********************************************************************************************************************/
	/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
	double Adash;


	if (launch == 2 || launch == 3) {
		//fp9 = fopen("Rinxy.csv", "w");

		if ((fp7 = fopen("FEM_LPmode_VCSELinput.csv", "r")) != NULL) {		//	VCSEL ��LP���[�h��1�������x���z�t�@�C�����J��

			int minput;
			char MPDfilename[256];
			char coupletype[128];
			if (couple == 0) { sprintf(coupletype, "V2F"); }
			else if (couple == 1) { sprintf(coupletype, "F2F"); }
			else { sprintf(coupletype, "error_warning"); }
			if (couple == 0) { max = num_mode_vin; }
			else if (couple == 1) { max = NLP; }

			for (int ii = 0; ii <= OFFrange / OFFres + 1; ii++) {
				if (launch == 2) { printf("\ncalculating MPD for lateral offset...\n"); }
				if (launch == 3) { printf("\ncalculating MPD for lateral offset of %d um...\n\n", ii * OFFres); }
				fscanf(fp7, "%s", trash);
				double E2m_evin, E2m_ev;
				for (minput = 0; minput < max; minput++) {
					if (ii == 0) {
						fscanf(fp7, "%d, %d, %lf, %lf, %lf, %lf, %lf,", &min, &lin, &tauin, &betain, &betain_devided_by_k, &Rinfin, &eigin);
						modemin[minput] = min;
						modelin[minput] = lin;
						Mbetain[minput] = betain;			}

					printf("%d\t%d\n", modemin[minput], modelin[minput]);
					if (couple == 0) { Adash = Avin; }		if (couple == 1) { Adash = A; }
					
					for (j = 0; j <= (int)(Adash / dr); j++) {
						if (ii == 0) {
							fscanf(fp7, "%lf,", &Rinlp[minput][j]);		}
						//printf("%lf, ", Rinlp[minput][j]);
					}
					if (ii == 0) { fscanf(fp7, "%lf\n", &Rinlp[minput][j]); }
					//printf("%lf\n", Rinlp[minput][j]);


					for (m = 0; m < NLP; m++) {

						Am_evev = 0.0;		Am_evod = 0.0;
						Am_odev = 0.0;		Am_odod = 0.0;
						E2m_evin = 0.0;		E2m_ev = 0.0;
						Rinxy = 0.0;

						for (i = 0; i < Nxy; i++) {
							xx1 = (-(double)Nxy*dx / 2.0) + ((double)i*dx);
							if (launch == 2) {
								xx2 = (r0 - (double)Nxy*dx / 2.0) + ((double)i*dx);
							}
							if (launch == 3) {
								xx2 = ((ii * (double)OFFres * 1.0e-6) - (double)Nxy*dx / 2.0) + ((double)i*dx);
							}

							for (j = 0; j < Nxy; j++) {
								yy = (-(double)Nxy*dx / 2.0) + ((double)j*dy);
								rr1 = sqrt(xx1*xx1 + yy * yy);
								rr2 = sqrt(xx2*xx2 + yy * yy);
								nrr1 = (int)(rr1 / dr);
								nrr2 = (int)(rr2 / dr);
								//printf("%lf\t%d\n", Rinxy, modemin[minput]);

								/*****   �����͓��˂��郌�[�U�i�t�@�C�o�j�́Cxy���ʏ��2�������x���z���v�Z    *****/
								if (nrr1 == 0) {		// �ȉ��̔z��Rinxy�́C���������̔z��Rlp���g���̂œǂݍ��݂�HDD�����葬��
									Rinxy = Rinlp[minput][0] + (Rinlp[minput][1] - Rinlp[minput][0]) * (rr1 / dr);
									//printf("%lf\t%d\n", Rinxy, modemin[minput]);
								}
								else if (nrr1 <= Nvin) {
									Rinxy = Rinlp[minput][nrr1] + (Rinlp[minput][nrr1 + 1] - Rinlp[minput][nrr1]) * ((rr1 - (double)nrr1*dr) / dr);
									//printf("%lf\t%d\n", Rinxy, modemin[minput]);
								}
								else if (nrr1 > Nvin) {
									//printf("%d, %lf, %d, %lf\n", minput, Mbetain[minput], Nvin, Avin);
									win = aa * sqrt(Mbetain[minput] * Mbetain[minput] - k * k*n0*n0);
									Rinxy = Rinlp[minput][Nvin] * (bessk(modemin[minput], win*rr1) / bessk(modemin[minput], win*Avin));
									//printf("%e, %e, %e, %e, %e,\t", win, Rinlp[minput][Nvin], bessk(minput, win*rr1)/bessk(minput, win*Avin), rr1, Avin);
									//printf("%lf, %d\n", Rinxy, nrr1);		//printf("%lf\t%d\n", Rinxy, modemin[minput]);
								}
								/*****   �����͎󂯎葤�̃t�@�C�o�ɂ�����Cxy���ʏ��2�������x���z���v�Z      *****/					/// �ȉ��̔z��Rinxy�́C���������̔z��Rlp���g���̂œǂݍ��݂�HDD�����葬��
								if (nrr2 == 0) {
									Rxy = Rlp[m][0] + (Rlp[m][1] - Rlp[m][0]) * (rr2 / dr);
								}
								else if (nrr2 <= N) {
									Rxy = Rlp[m][nrr2] + (Rlp[m][nrr2 + 1] - Rlp[m][nrr2]) * ((rr2 - (double)nrr2*dr) / dr);
								}
								else if (nrr2 > N) {
									w = aa * sqrt(Mbeta[m] * Mbeta[m] - k * k*n0*n0);
									Rxy = Rlp[m][N] * (bessk(modem[m], w*rr2) / bessk(modem[m], w*A));
								}

								/*****   �d�Ȃ�ϕ��ɗp����d�ꕪ�z�̎Z�o�E�d�Ȃ�ϕ��̎��s   *****/
								/*if (minput == 1) {
									if (m == 1) {
										printf("%lf, %lf\t%d\n", Rinxy, Rxy, modemin[minput]); }}*/
										//fprintf(fp9, "%lf,", Rinxy);
										/*if (Rinxy == Rxy) { printf("Y "); }
										else printf("N ");*/
								Em_evin = Rinxy * cos((double)(modemin[minput]) * atan2(yy, xx1));
								Em_odin = Rinxy * sin((double)(modemin[minput]) * atan2(yy, xx1));
								Em_ev = Rxy * cos((double)(modem[m]) * atan2(yy, xx2));
								Em_od = Rxy * sin((double)(modem[m]) * atan2(yy, xx2));
								if (i < 1 && j < 1) {
									/*if (modemin[minput] == modem[m]) {
										//printf("%lf, %lf, %lf, %lf\t\t", Em_evin, Em_odin, Em_ev, Em_od); 
										printf("%d   ", m);
										printf("%d, %d, %lf \t %d, %d, %lf\n", modemin[minput], modelin[minput], Mbetain[minput], modem[m], model[m], Mbeta[m]);
									}*/
									if (modem[m] == 0) {
										//printf("%lf, %lf, %lf, %lf\t\t", Em_evin, Em_odin, Em_ev, Em_od); 
										printf("%d ", m);
										//printf("%d, %d, %lf \t %d, %d, %lf\n", modemin[minput], modelin[minput], Mbetain[minput], modem[m], model[m], Mbeta[m]);
									}
									else { printf("%d ", m); }
								}

								E2m_evin = E2m_evin + Em_evin * Em_evin*dx*dy;
								E2m_ev = E2m_ev + Em_ev * Em_ev*dx*dy;
								Am_evev = Am_evev + Em_evin * Em_ev*dx*dy;
								Am_evod = Am_evod + Em_evin * Em_od*dx*dy;
								Am_odev = Am_odev + Em_odin * Em_ev*dx*dy;
								Am_odod = Am_odod + Em_odin * Em_od*dx*dy;
							}

						}

						//printf("%lf\t", E2m_evin);

						cef_evev = Am_evev * Am_evev / (E2m_evin * E2m_ev);
						cef_evod = Am_evod * Am_evod / (E2m_evin * E2m_ev);
						cef_odev = Am_odev * Am_odev / (E2m_evin * E2m_ev);
						cef_odod = Am_odod * Am_odod / (E2m_evin * E2m_ev);

						/*Am_evev = Am_evev * Am_evev / (((2.0*omega*Mu0) / Mbetain[minput]) * (2.0*omega*Mu0) / Mbeta[m]);
						Am_evod = Am_evod * Am_evod / (((2.0*omega*Mu0) / Mbetain[minput]) * (2.0*omega*Mu0) / Mbeta[m]);
						Am_odev = Am_odev * Am_odev / (((2.0*omega*Mu0) / Mbetain[minput]) * (2.0*omega*Mu0) / Mbeta[m]);
						Am_odod = Am_odod * Am_odod / (((2.0*omega*Mu0) / Mbetain[minput]) * (2.0*omega*Mu0) / Mbeta[m]);*/
						if (modem[m] == 0 && modemin[minput] == 0) {
							MPD2d[m][minput] = 100.0*(cef_evev + cef_evod + cef_odev + cef_odod);
						}
						if (modem[m] != 0 && modemin[minput] == 0) {
							MPD2d[m][minput] = 100.0*(cef_evev + cef_evod + cef_odev + cef_odod) / 2;
						}
						if (modem[m] == 0 && modemin[minput] != 0) {
							MPD2d[m][minput] = 100.0*(cef_evev + cef_evod + cef_odev + cef_odod) / 2;
						}
						if (modem[m] != 0 && modemin[minput] != 0) {
							MPD2d[m][minput] = 100.0*(cef_evev + cef_evod + cef_odev + cef_odod) / 2;
						}
						/*if (modem[m] == 0 & modemin[minput]==0) {
							MPD2d[m][minput] = 100.0*(Am_evev + Am_evod + Am_odev + Am_odod);	}
						if (modem[m] != 0 & modemin[minput] == 0){
							MPD2d[m][minput] = 100.0*(Am_evev + Am_evod + Am_odev + Am_odod) / 2;	}
						if (modem[m] == 0 & modemin[minput] != 0){
							MPD2d[m][minput] = 100.0*(Am_evev + Am_evod + Am_odev + Am_odod) / 2;	}
						if (modem[m] != 0 & modemin[minput] != 0) {
							MPD2d[m][minput] = 100.0*(Am_evev + Am_evod + Am_odev + Am_odod) / 4;	}*/

							//printf("%lf, %lf, %lf, %lf\n", Am_evev*100, Am_evod*100, Am_odev*100, Am_odod*100);
					}
					printf("\n");
				}


				printf("\n");
				for (minput = 0; minput < max; minput++) { printf("%d\t%lf\n", modemin[minput], Mbetain[minput]); }
				printf("\n");
				for (m = 0; m < NLP; m++) { printf("%d\t%lf\n", modem[m], Mbeta[m]); }

				if (launch == 2) { sprintf(MPDfilename, "MPD_%s\\MPD2d_%s_%dum_Offset_single(not_dependence_on_x).csv", coupletype, coupletype, r0); }
				if (launch == 3) { sprintf(MPDfilename, "MPD_%s\\MPD2d_%s_%dum_Offset.csv", coupletype, coupletype, ii*OFFres); }


				if ((fp8 = fopen(MPDfilename, "w")) != NULL) {
					fprintf(fp8, "Mode number for receiver-side fiber,");
					if (couple == 0) { fprintf(fp8, "m, l, "); }
					for (j = 0; j < max; j++) { fprintf(fp8, "%d,", j); }
					fprintf(fp8, "\n");
					for (i = 0; i < NLP; i++) {
						fprintf(fp8, "%d,", i);
						if(couple == 0) { fprintf(fp8, "%d, %d, ", modem[i], model[i]); }
						for (j = 0; j < max; j++) {
							// printf("%e, ", MPD2d[i][j]);
							fprintf(fp8, "%e,", MPD2d[i][j]);
						}
						fprintf(fp8, "\n");
					}
					fclose(fp8);
				}
				if (launch == 2) break; 
			}
		}
	}	// if(launch==2) �̏I��
		/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
		/*********************************************************************************************************************/
	for (m = 0;m < 80;m++) {
		printf("%d, ", modem[m]);
	}



	printf("\n���[�h�p���[���z�v�Z����\n");

	/** 5. �C���p���X�����̌v�Z **/
	taumax = taumin = Mtau[0];
	for (myu = 0; myu < NLP; myu++) {
		if (taumax < Mtau[myu]) { taumax = Mtau[myu]; }
		if (taumin > Mtau[myu]) { taumin = Mtau[myu]; }
	}

	Ntmax = int((taumax*L) / Tv), Ntmin = int((taumin*L) / Tv);
	P = drealvector(0, Ntmax - Ntmin); init_realvector(P, 0, Ntmax - Ntmin); //�x���X�e�b�v���Ƃ̃p���[�i�[�z��
	for (myu = 0; myu < NLP; myu++) {
		x = int((Mtau[myu] * L) / Tv) - Ntmin;
		//P[x] = P[x] +1.0;}
		if (launch == 0) {
			if (modem[myu] == 0) { P[x] = P[x] + 2.0; }
			else { P[x] = P[x] + 4.0; }
		}
		if (launch == 1) {
			//if ( modem[myu] == 1 ) { P[x] = P[x] + 2.0*MPD[myu]; } else { P[x] = P[x] + 4.0*MPD[myu]; }}
			P[x] = P[x] + MPD[myu];
		}
	}

	for (x = 0; x <= Ntmax - Ntmin; x++) { fprintf(fp5, "%lf,%lf\n", double(Ntmin + x)*Tv, P[x]); }
	printf("�C���p���X�����v�Z����\n");

	/** 6. �`�B�֐�Hw�̌v�Z **/
	/* ���g������ M */
	for (y = 0; y <= Nf; y++) {
		Hw = ReHw = ImHw = bw = 0.0;
		for (x = 0; x <= Ntmax - Ntmin; x++) {
			ReHw = ReHw + P[x] * cos((2.0*PI*(fmin + (double)y*dfrq))*(double)x*Tv);
			ImHw = ImHw - P[x] * sin((2.0*PI*(fmin + (double)y*dfrq))*(double)x*Tv);
		}

		M[y] = sqrt(ReHw*ReHw + ImHw * ImHw) / double(Ntotal);
		fprintf(fp6, "%lf,", M[y]);
	}

	/** -3dB�ш� bw **/
	if (M[Nf] < 0.5) {
		for (y = 0; y <= Nf - 1; y++) {
			if (M[y] > 0.5 && M[y + 1] < 0.5) { bw = (fmin*1.0e3) + (dfrq *1.0e3)*(y + (M[y] - 0.5) / (M[y] - M[y + 1])); break; }
		}
	}
	fprintf(fp6, "%lf,", bw);

	printf("���g�������v�Z����\n");

	/** 7. �o�˔g�`�̌v�Z **/
	LT = 2000;
	xmax = Ntmax - Ntmin + 1;
	xL = sizeof(double)*(xmax + LT);
	opulse = drealvector(0, xL); init_realvector(opulse, 0, xL);
	opulse = convolution(xmax, P, LT, pulse);
	for (i = 0; i <= xmax + LT + 1; i++) { fprintf(fopulse, "%lf\n", opulse[i]); }

	/* �z��̋L���̈��� */
	free_drealvector(GI, 0);
	free_drealvector(q, 0);
	free_drealvector(qg, 0);
	free_drealvector(R, 0);
	free_drealvector(R2, 0);
	free_drealvector(Rb, 0);
	free_drealvector(a, 0);
	free_drealvector(b, 0);
	free_drealvector(ML, 0);
	free_drealvector(MD, 0);
	free_drealvector(Mtau, 0);
	free_dintvector(modem, 0);
	free_drealvector(P, 0);
	free_drealvector(M, 0);
	//
	free_drealvector(Mbeta, 0);
	free_drealvector(MPD, 0);
	free_drealvector(pulse, 0);
	free_drealvector(opulse, 0);

	fprintf(fp2, "Total number of LP modes (WKB), Nwkb,%d\n", Nwkb);
	fprintf(fp2, "Total number of LP modes, NLP,%d\n", NLP);
	fprintf(fp2, "Total number of LP modes with 0th m-order, NLP0,%d\n", NLP0);
	fprintf(fp2, "Total number of modes, Ntotal,%d\n", Ntotal);

	fclose(fp2);
	fclose(fp4);
	fclose(fp5);
	fclose(fp6);	//fclose(fp7);
	fclose(fopulse);
	fclose(fmpd);
	//free_dmatrix(Rlp, 0, sizeofRLP, 0, sizeofRLP);	//free_dmatrix(Rinlp, 0, (int)(2 * A / dr), 0, (int)(2 * A / dr));
	//free_dmatrix(MPD2d, 0, NLP, 0, NLP);

	system("pause");

	return 0;
}



/*A.1. ��2��ό`�x�b�Z���֐� bessk (n, x) */
// ��1��ό`Bessel�֐��in=0�jI0(x)
double bessi0(double x)
{
	double ax, ans, y;
	// Polynomial fit
	if ((ax = fabs(x)) < 3.75) {
		y = x / 3.75;
		y *= y;
		ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
			+ y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
	}
	else {
		y = 3.75 / ax;
		ans = (exp(ax) / sqrt(ax))*(0.39894228 + y * (0.1328592e-1
			+ y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2
				+ y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1
					+ y * 0.392377e-2))))))));
	}
	return ans;
}

// ��2��ό`Bessel�֐��in=0�j K0(x)
double bessk0(double x)
{
	double bessi0(double x);
	double y, ans;
	//polynomial fit
	if (x <= 2) {
		y = x * x / 4.0;
		ans = (-log(x / 2.0)*bessi0(x)) + (-0.57721566 + y * (0.42278420
			+ y * (0.23069756 + y * (0.3488590e-1 + y * (0.262698e-2
				+ y * (0.10750e-3 + y * 0.74e-5))))));
	}
	else {
		y = 2.0 / x;
		ans = (exp(-x) / sqrt(x))*(1.25331414 + y * (-0.7832358e-1
			+ y * (0.2189568e-1 + y * (-0.1062446e-1 + y * (0.587872e-2
				+ y * (-0.251540e-2 + y * 0.53208e-3))))));
	}
	return ans;
}

// ��1��ό`Bessel�֐��in=1�j I1(x)
double bessi1(double x)
{
	double ax, ans, y;
	if ((ax = fabs(x)) < 3.75) {
		y = x / 3.75;	   y *= y;
		ans = ax * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934
			+ y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
	}
	else {
		y = 3.75 / ax;
		ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1
			- y * 0.420059e-2));
		ans = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2
			+ y * (0.163801e-2 + y * (-0.1031555e-1 + y * ans))));
		ans *= (exp(ax) / sqrt(ax));
	}
	return x < 0.0 ? -ans : ans;
}

// ��2��ό`Bessel�֐��in=1�j K1 (x)
double bessk1(double x) {
	double bessi1(double x);
	double y, ans;
	if (x <= 2.0) {
		y = x * x / 4.0;
		ans = (log(x / 2.0)*bessi1(x)) + (1.0 / x)*(1.0 + y * (0.15443144
			+ y * (-0.67278579 + y * (-0.18156897 + y * (-0.1919402e-1
				+ y * (-0.110404e-2 + y * (-0.4686e-4)))))));
	}
	else {
		y = 2.0 / x;
		ans = (exp(-x) / sqrt(x))*(1.25331414 + y * (0.23498619
			+ y * (-0.3655620e-1 + y * (0.1504268e-1 + y * (-0.780353e-2
				+ y * (0.325614e-2 + y * (-0.68245e-3)))))));
	}
	return ans;
}

// ��2��ό`Bessel�֐� Kn(x)
double bessk(int n, double x) {
	double bessk0(double x);
	double bessk1(double x);
	double bk, bkm, bkp, tox;   int j;
	if (n == 0) return bessk0(x);
	if (n == 1) return bessk1(x);
	if (n >= 2) {
		tox = 2.0 / x;
		bkm = bessk0(x);
		bk = bessk1(x);
		for (j = 1; j < n; j++) {
			bkp = bkm + j * tox*bk;
			bkm = bk;
			bk = bkp;
		}
		return bk;
	}
	else return 0;
}

/*A.2. ���ܗ��̔g�������֐�dndl (n_lamda, mater) */
double dndl(double lamda, double n_lamda, int mater)
{
	if (mater == 0) // PMMA
		return ((0.02173e-4*lamda - 18.79107e-4)*n_lamda + (-0.03109e-4*lamda + 26.85035e-4))*1.0e9;
	else if (mater == 1) // CYTOP
		return ((0.02383e-5*lamda - 25.34898e-5)*n_lamda + (-0.02983e-5*lamda + 31.46375e-5))*1.0e9;
	else if (mater == 2) //Silica
		return ((0.02817e-5*lamda - 29.61218e-5)*n_lamda + (-0.03807e-5*lamda + 39.02595e-5))*1.0e9;
	else if (mater == 3) //PTCEMA
		return ((0.0076424e-4*lamda - 7.92635e-4)*n_lamda + (-0.0107534e-4*lamda + 11.10000e-4))*1.0e9;
	else printf(" Material number is not correct !\n"); exit(EXIT_FAILURE);					}



/*A.3. �W���s��v�f�v�Z�֐� S_matrix (a, b, q, m, n, v, w, D)*/
/********************************************************************************************************************************************/
/// b[0]=S00, b[1]=S11, ...................... b[j]=Sjj ......  b[n]=Snn					// a[0]=___, a[1]=S01, a[2]=S12, ... a[j]=Sj-1,j ... a[n]=Sn-1,n					// a[0]=___, a[1]=S10, a[2]=S21, ... a[j]=Sj,j-1 ... a[n]=Sn,n-1
void S_matrix(double *a, double *b, double *q, int m, int n, double v, double w, double D) {
	if (m == 0) {
		b[0] = -(1.0 / 2.0) + (3.0*q[0] + 2.0*q[1])*(v*v / 60.0)*((D*D) / (n*n)) - (1.0 / 12.0)*(w*w)*((D*D) / (n*n));
		b[n] = -((2.0*n - 1.0) / 2.0) + ((5.0*n - 2.0)*q[n - 1] + 3.0*(5.0*n - 1.0)*q[n])*(v*v / 60.0)*((D*D) / (n*n)) - ((4.0*n - 1.0) / 12.0)*(w*w)*((D*D) / (n*n))
			- w * D*bessk(1, w*D) / bessk(0, w*D);
	}
	else {
		b[0] = 0.0;
		b[n] = -(1.0 - m * m)*((2.0*n - 1.0) / 2.0) + ((5.0*n - 2.0)*q[n - 1] + 3.0*(5.0*n - 1.0)*q[n])*(v*v / 60.0)*((D*D) / (n*n)) - ((4.0*n - 1.0) / 12.0)*(w*w)*((D*D) / (n*n))
			- (m*m)*((n - 1.0)*(n - 1.0))*log((double)n / ((double)n - 1.0)) - m * m - w * D*bessk(m - 1, w*D) / bessk(m, w*D) - m;
	}
	///----------------------------------------
	b[1] = -2.0*(1.0 - m * m) + (3.0*q[0] + 30 * q[1] + 7.0*q[2])*(v*v / 60.0)*((D*D) / (n*n)) - (2.0 / 3.0)*(w*w)*((D*D) / (n*n)) - (m*m)*4.0*log(2.0);			int j;
	///----------------------------------------
	for (j = 2; j < n; j++) {
		b[j] = -2.0*j*(1.0 - m * m) + ((5.0*j - 2.0)*q[j - 1] + 30.0*j*q[j] + (5.0*j + 2.0)*q[j + 1])*(v*v / 60.0)*((D*D) / (n*n)) - (2.0 / 3.0)*j*(w*w)*((D*D) / (n*n))
			- (m*m)*(((j - 1.0)*(j - 1.0))*log((double)j / ((double)j - 1.0)) + ((j + 1.0)*(j + 1.0))*log(((double)j + 1.0) / (double)j));
	}
	///----------------------------------------
	a[0] = 0.0;//���g�p�v�f�ɂ�0���i�[
	a[1] = (1.0 / 2.0) + (2 * q[0] + 3 * q[1])*(v*v / 60.0)*((D*D) / (n*n)) - (1.0 / 12.0)*(w*w)*((D*D) / (n*n)) - (m*m) / 2;
	///----------------------------------------
	for (j = 2; j <= n; j++) {
		a[j] = ((2.0*j - 1.0) / 2.0)*(1.0 - m * m) + ((5.0*j - 3.0)*q[j - 1] + (5.0*j - 2.0)*q[j])*(v*v / 60.0)*((D*D) / (n*n)) - ((2.0*j - 1.0) / 12.0)*(w*w)*((D*D) / (n*n))
			+ (m*m)*(j - 1.0)*j*log((double)j / ((double)j - 1.0));
	}
}
/********************************************************************************************************************************************/


/*A.4. �����R���X�L�[����	mcholesky ( a, b, ML, MD, m, n )*/
void mcholesky(double *a, double *b, double *ML, double *MD, int m, int n) {
	int i;
	if (m == 0) {
		ML[0] = 0.0;//���g�p�v�f�ɂ�0���i�[
		MD[0] = b[0];
		for (i = 1; i <= n; i++) {
			MD[i] = b[i] - a[i] * a[i] / MD[i - 1];
			ML[i] = a[i] / MD[i - 1];
		}
	}
	else {
		ML[0] = 0.0;//���g�p�v�f�ɂ�0���i�[
		ML[1] = 0.0;
		MD[0] = 0.0;
		MD[1] = b[1];
		for (i = 2; i <= n; i++) {
			MD[i] = b[i] - a[i] * a[i] / MD[i - 1];
			ML[i] = a[i] / MD[i - 1];
		}
	}
}


/*A.5. �����R���X�L�[����@�ɂ�鋁��	mcholesky_sol ( a, b, ML, MD, m, n )*/
void mcholesky_sol(double *ML, double *MD, double *R, int m, int n) {
	int i;
	//�uLy=R0�v������
	if (m == 0) {
		for (i = 1; i <= n; i++) {
			R[i] = R[i] - R[i - 1] * ML[i];
		}
		//�u(D(LT))R1=y�v���u(LT)R1=(D-1)y=y'�v�ɕς���
		for (i = 1; i <= n; i++) {
			R[i] = R[i] / MD[i];
		}
		//�u (LT)R1=y'�v������
		for (i = n - 1; i >= 0; i--) {
			R[i] = R[i] - ML[i + 1] * R[i + 1];
		}
	}

	else {
		for (i = 2; i <= n; i++) {
			R[i] = R[i] - R[i - 1] * ML[i];
		}
		//�u(D(LT))R1=y�v���u(LT)R1=(D-1)y=y'�v�ɕς���
		for (i = 2; i <= n; i++) {
			R[i] = R[i] / MD[i];
		}
		//�u (LT)R1=y'�v������
		for (i = n - 1; i >= 1; i--) {
			R[i] = R[i] - ML[i + 1] * R[i + 1];
		}
	}
}

/*A.6. �t�ׂ���@�̏����x�N�g���v�Z	R0 ( MD, R, m, n )*/
/* �Ίp�s��D�̐������ő�ƂȂ�v�f����1�ł���悤�ȃx�N�g����I��*/
void R0(double *MD, double *R, int m, int n) {
	int i, j = 1;
	if (m == 0) {
		for (i = 0; i <= n - 1; i++) {
			if (fabs(MD[i + 1]) < fabs(MD[j])) { j = i + 1; }
		}
	}
	else {
		for (i = 1; i <= n - 1; i++) {
			if (fabs(MD[i + 1]) < fabs(MD[j])) { j = i + 1; }
		}
	}
	//R0�̏����l�̑��
	for (i = 0; i <= n; i++) {
		if (i == j) { R[i] = 1.0; }
		else { R[i] = 0.0; }
	}
}

/*A.7. �t�ׂ���@�̉��x�N�g���K�i��	R_norm ( R, n )*/
void R_norm(double *R, int n) {
	int i;	double s = 0;
	for (i = 0; i <= n; i++) { s = s + R[i] * R[i]; }			//�s��v�f�̂Q��a
	if (s != 0) {
		for (i = 0; i <= n; i++) { R[i] = R[i] / sqrt(s); }
	}
}

/*A.8. �ŗL�l�v�Z�iRayleigh quotient�jEigen ( R, a, b, m, n )*/// R�x�N�g�����K�i�����Ă��邽�ߓ��ς�1
double Eigen(double *R, double *a,
	double *b, int n, int m) {
	int i;	double s = 0;
	if (m == 0) {
		s = (R[0] * b[0] + R[1] * a[1])*R[0];
		for (i = 1; i < n; i++) {
			s += (R[i - 1] * a[i] + R[i] * b[i] + R[i + 1] * a[i + 1])*R[i];
		}
		s += (R[n - 1] * a[n] + R[n] * b[n])*R[n];
		return s;
	}
	else {
		s = (R[1] * b[1] + R[2] * a[2])*R[1];
		for (i = 2; i < n; i++) {
			s += (R[i - 1] * a[i] + R[i] * b[i] + R[i + 1] * a[i + 1])*R[i];
		}
		s += (R[n - 1] * a[n] + R[n] * b[n])*R[n];
		return s;
	}
}

/*A.9. �Q�x���v�Z�p�֐�	dbdk_bunbo ( R, D, w, m, n )*/
/* ���̓p�����[�^�i�������d�ꕪ�z�C�R�A�a�C�`���萔�C�v�f�������j�ɑ΂���m�����[�h�Q�x���v�Z���̕��ꕪ�q���v�Z���� */

// �������d�ꐬ�� R[0]~R[n], �K�i���`���萔 w, �K�i���R�A�a D, ���ʊp���[�h���� m, ������ n
double dbdk_bunbo(double *R, double D, double w, int m, int n) {
	int i;	double s = 0;
	for (i = 0; i <= n - 1; i++) {
		s = s + (1.0 / 12.0) * ((D*D) / (n*n)) * ((double)(4 * i + 1)*R[i] * R[i] + 2.0*(double)(2 * i + 1)*R[i] * R[i + 1] + (double)(4 * i + 3)*R[i + 1] * R[i + 1]);
	}
	if (m == 0) {
		return s + ((bessk(1, w*D)*bessk(1, w*D) / (bessk(0, w*D)*bessk(0, w*D))) - 1.0) * ((D*D)*(R[n] * R[n]) / 2.0);
	}
	else {
		return s + ((bessk(m - 1, w*D)*bessk(m + 1, w*D) / (bessk(m, w*D)*bessk(m, w*D))) - 1.0) * ((D*D)*(R[n] * R[n]) / 2.0);
	}
}

/*A.10. �Q�x���v�Z�p�֐�	dbdk_bunshi ( R,qg, D, w, m, n )*/
/* ���̓p�����[�^�i�������d�ꕪ�z�C�R�A�a�C�`���萔�C�v�f�������j�ɑ΂���m�����[�h�Q�x���v�Z���̕��ꕪ�q���v�Z���� */

// �������d�ꐬ�� R[0]~R[n], �K�i���`���萔 w, �K�i���R�A�a D, ���ʊp���[�h���� m, ������ n
// ���ܗ����U�p�����[�^ qg[0]~qg[n] ( = n*(d(kn)/dk) )
double dbdk_bunshi(double *R, double *qg, double D, double w, int m, int n)
{
	int i;		double s = 0;
	for (i = 0; i <= n - 1; i++) {
		s = s + (1.0 / 12.0) * ((D*D) / (n*n))* (((double)(3 * i) + 3.0 / 5.0)*qg[i] * R[i] * R[i] + ((double)i + 2.0 / 5.0)*(2.0*qg[i] * R[i + 1] + qg[i + 1] * R[i])*R[i] + ((double)i + 3.0 / 5.0)*(qg[i] * R[i + 1] + 2.0*qg[i + 1] * R[i])*R[i + 1] + ((double)(3 * i) + 12.0 / 5.0)*qg[i + 1] * R[i + 1] * R[i + 1]);
	}
	if (m == 0) {
		return s + qg[n] * ((bessk(1, w*D)*bessk(1, w*D) / (bessk(0, w*D)*bessk(0, w*D))) - 1.0) * ((D*D)*(R[n] * R[n]) / 2.0);
	}
	else {
		return s + qg[n] * ((bessk(m - 1, w*D)*bessk(m + 1, w*D) / (bessk(m, w*D)*bessk(m, w*D))) - 1.0) * ((D*D)*(R[n] * R[n]) / 2.0);
	}
}

/* A.11. ��ݍ��ݐϕ��֐� convolution ( xmax, P, LT, pulse) */
double *convolution(int xmax, double *P, int LT, double *pulse) {
	int i, j, xL;
	double *R;
	xL = sizeof(double)*(xmax + LT + 1);
	R = (double*)malloc(xL);
	printf("xmax = %d\n" "LT = %d\n" "xL = %d\n", xmax, LT, xL);
	//for ( i = 0 ; i <= xmax ; i++ ){ printf ( "%f\n", P[i] );}
	//for ( i = 0 ; i <= LT ; i++ ){ printf ( "%f\n", pulse[i] );}
	for (j = 0; j <= xmax + LT; j++) { R[j] = 0; }
	for (i = 0; i < xmax; i++) {
		for (j = 0; j <= LT; j++) {
			R[i + j] = R[i + j] + P[i] * pulse[j];
			//printf ( "%f\n", R[i+j] );
		}
	}				//for ( i = 0; i <= xmax+LT+1 ; i++ ){ printf ( "%f\n", R[i] );}
	return R;
}

/* A.12. �����x�N�g���̈�m�ۗp�֐� dvector (i, j) */
int *dintvector(int i, int j) {
	int *a;
	if ((a = (int *)malloc((j - i + 1) * sizeof(int))) == NULL) {
		printf("Memory cannot be allocated !\n"); exit(1);
	}
	return (a - i);
}

/* A.13. �����x�N�g���̈����p�֐� free_dvector (a, i) */
void free_dintvector(int *a, int i) {
	free((void*)(a + i));
}


/* A.14. �����x�N�g���̈�m�ۗp�֐� dvector (i, j) */
double *drealvector(int i, int j) {
	double *a;
	if ((a = (double *)malloc((j - i + 1) * sizeof(double))) == NULL) {
		printf("Memory cannot be allocated !\n"); exit(1);
	}
	return (a - i);
}

/* A.15. �����x�N�g���̈����p�֐� free_dvector (a, i) */
void free_drealvector(double *a, int i) {
	free((void*)(a + i));
}


/* A.16. �����s��̈�m�ۗp�֐� dmatrix ( nr1, nr2, nl1, nl2 ) */
double **dmatrix(int nr1, int nr2, int nl1, int nl2) {
	// nrow: �s�̐�, ncol: ��̐�
	double **a;
	int i, nrow, ncol;
	nrow = nr2 - nr1 + 1;
	ncol = nl2 - nl1 + 1;
	/* �s�̊m�� */
	if ((a = (double **)malloc(nrow * sizeof(double*))) == NULL) {
		printf("Memory cannot be allocated !\n"); exit(1);
	}
	a = a - nr1; // �s�����炷
				 /* ��̊m�� */
	for (i = nr1; i <= nr2; i++) a[i] = (double *)malloc(ncol * sizeof(double));
	for (i = nr1; i <= nr2; i++) a[i] = a[i] - nl1;	// ������炷
	return (a);
}

/* A.17.�����s��̈����p�֐� free_dmatrix ( a, nr1, nr2, nl1, nl2 ) */
void free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2) {
	int i;
	for (i = nr1; i <= nr2; i++) free((void*)(a[i] + nl1));
	free((void*)(a + nr1));
}


/* A.18. �����x�N�g���������֐� init_vector ( a, nr1, nr2 ) */
void init_intvector(int *a, int nr1, int nr2) {
	int i;
	for (i = nr1; i <= nr2; i++) { a[i] = 0; }
}


/* A.19. �����x�N�g���������֐� init_vector ( a, nr1, nr2 ) */
void init_realvector(double *a, int nr1, int nr2) {
	int i;
	for (i = nr1; i <= nr2; i++) { a[i] = 0.0; }
}


/* A.21. �G���[�֐� -------(A.22��A.23�ɗp����)-----------*/
void nrerror(char error_text[]) {
	fprintf(stderr, "Numerical Recipes run=time error///\n");
	fprintf(stderr, "%s\n", error_text);
	fprintf(stderr, "...now exiting to system...\n");
	exit(1);
}


/* A.22. float�^��2�����z�� matrix(long nrl, long nrh, long ncl, long nch) */
float **matrix(long nrl, long nrh, long ncl, long nch) {
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	float **m;
	m = (float **)malloc((size_t)((nrow + NR_END) * sizeof(float*)));
	// if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;		  m -= nrl;
	m[nrl] = (float *)malloc((size_t)((nrow*ncol + NR_END) * sizeof(float*)));
	// if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;	  m[nrl] -= ncl;
	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;
	return m;
}

/* A.23. float�^2�����z��̉�� matrix(long nrl, long nrh, long ncl, long nch) */
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch) {
	free((FREE_ARG)(m[nrl] + ncl - NR_END));
	free((FREE_ARG)(m + nrl - NR_END));
}

/* A.24. float�^2�����z��̉�� f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh) */
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
	long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
	double ***t;
	/* allocate pointers to pointers to rows */
	t = (double ***)malloc((size_t)((nrow + NR_END) * sizeof(double**)));
	//if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;
	/* allocate pointers to rows and set pointers to them */
	t[nrl] = (double **)malloc((size_t)((nrow*ncol + NR_END) * sizeof(double*)));
	//if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	/* allocate rows and set pointers to them */
	t[nrl][ncl] = (double *)malloc((size_t)((nrow*ncol*ndep + NR_END) * sizeof(double)));
	//if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;
	for (j = ncl + 1;j <= nch;j++) t[nrl][j] = t[nrl][j - 1] + ndep;
	for (i = nrl + 1;i <= nrh;i++) {
		t[i] = t[i - 1] + ncol;
		t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
		for (j = ncl + 1;j <= nch;j++) t[i][j] = t[i][j - 1] + ndep;
	}
	/* return pointer to array of pointers to rows */
	return t;
}

/* A.24. float�^2�����z��̉�� f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh) */
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
	/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
	free((FREE_ARG)(t[nrl] + ncl - NR_END));
	free((FREE_ARG)(t + nrl - NR_END));
}

/* A.25. �����s�񏉊����֐� init_dmatrix ( a, nr1, nr2, nl1, nl2 ) */
void init_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2) {
	for (int i = nr1; i <= nr2; i++) {
		for (int j = nl1; j <= nl2; j++) { a[i][j] = 0.0; }
	}
}

/* A.26. �����s�񏉊����֐� init_d3tensor ( a, nr1, nr2, nl1, nl2, np1, np2 ) */
void init_d3tensor(double ***a, int nr1, int nr2, int nl1, int nl2, int np1, int np2) {
	for (int i1 = nr1; i1 <= nr2; i1++) {
		for (int i2 = nl1; i2 <= nl2; i2++) { 
			for(int i3 = np1; i3 <= np2; i3++) { a[i1][i2][i3] = 0.0; }
}}}

/* A.27. �����s�񏉊����֐� init_matrix ( a, nr1, nr2, nl1, nl2 ) */
void init_matrix(float **a, int nr1, int nr2, int nl1, int nl2) {
	for (int i = nr1; i <= nr2; i++) {
		for (int j = nl1; j <= nl2; j++) { a[i][j] = 0.0; }
	}
}