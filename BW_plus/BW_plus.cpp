#pragma warning( disable: 0266 )	//!? �ufprintf�������܂��ł��v�����
#pragma warning( disable: 4101 )	//!? �g�p���Ă��Ȃ��ϐ�
#pragma warning( disable: 4477 )	//!? ����������Ă��Ȃ��ϐ��ɑ΂���fscanf�Œl����
#pragma warning( disable: 4700 )	//!? ����������Ă��Ȃ��z��ɑ΂���fscanf�Œl����
#pragma warning( disable: 4996 )	//!? fscanf�֘A���āA�o�b�t�@�I�[�o�[�������N�����Ȃǂ̐Ǝ㐫�̋���
#pragma warning( disable: 6001 )	//!? ����������Ă��Ȃ����������g�p����Ă��܂�
#pragma warning( disable: 6031 )	//!? fscanf�̖߂�l�𖳎�
#pragma warning( disable: 6054 )	//!? strcmp�֐��Ɋւ���x���𖳎�
#pragma warning( disable: 6385 )	
#pragma warning( disable: 6386 )	
#pragma warning( disable: 6262 )	
#pragma warning( disable: 26451 )	

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <direct.h>
#include <Windows.h>
#include <tchar.h>
#include <locale>
#include <cstdlib>
#define _USE_MATH_DEFINES
#define NR_END 1
#define FREE_ARG char*
using namespace std;

//! ���萔�̒�`
#define PI 3.141592653
#define C 2.99792458e8
#define Epsilon0 8.8541878e-12
#define Mu0 PI*4.0e-7

//! �֐����v���O�����̓ǂݍ���
// ���ܗ��g������
double dndl ( double lamda, double n_lamda, int mater );
// ���ܗ��Z�x����
double dndc ( double lamda, int mater );
// �R�A���S���ܗ�
double ncore ( double lamda, int mater );
// �N���b�h���ܗ�
double nclad ( double lamda, int mater );
// ��1��ό`Bessel�֐� In(x)
double  bessi0 (double x), bessi1 (double x);
// ��2��ό`Bessel�֐� Kn(x)
double  bessk0 (double x),  bessk1 (double x), bessk ( int n, double x );
// �W���s��v�f�v�Z�֐�
void S_matrix ( double *a, double *b, double *q, int m, int n, double v, double w, double D );
// �����R���X�L�[����
void mcholesky ( double *a, double *b, double *ML, double *MD, int m, int n );
// �����R���X�L�[����@�ɂ�������������
void mcholesky_sol ( double *ML, double *MD, double *R, int m, int n );
// �t�ׂ���@�̏����x�N�g���v�Z
void R0 ( double *MD, double *R, int m, int n );
// �t�ׂ���@�̉��x�N�g���K�i��
void R_norm ( double *R, int n );
// �ŗL�l�v�Z
double Eigen ( double *R, double *a, double *b, int n, int m );
// �Q�x���v�Z�p�֐�
double dbdk_bunbo ( double *R, double D, double w, int m, int n );
// �Q�x���v�Z�p�֐�
double dbdk_bunshi ( double *R, double *q2, double D, double w, int m, int n);
// �������̈�m�ہi�����x�N�g���p�j
int *dintvector ( int i, int j );
// �������̈����i�����x�N�g���p�j
void free_dintvector ( int *a,  int i );
// �������̈�m�ہi�����x�N�g���p�j
double *drealvector ( int i, int j );
// �������̈����i�����x�N�g���p�j
void free_drealvector ( double *a,  int i );
// �������̈�m�ہi�}�g���N�X�p�j
double **dmatrix ( int nr1, int nr2, int nl1, int nl2 ); 
// �������̈����i�}�g���N�X�p�j
void free_dmatrix ( double **a, int nr1, int nr2, int nl1, int nl2 );
// 2���������z�񏉊���
void init_dmatrix ( double **a, int nr1, int nr2, int nl1, int nl2 );
// 1���������z�񏉊���
void init_intvector ( int *a, int nr1, int nr2 );
// 1���������z�񏉊���
void init_realvector ( double *a, int nr1, int nr2 );
// 3�����z��
double*** d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void init_d3tensor(double*** a, int nr1, int nr2, int nl1, int nl2, int np1, int np2);
void free_d3tensor(double*** t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
// ��ݍ��ݐϕ����s�֐�
double* convolution ( int n1, double* P1, int n2, double* P2 );
// �f�B���N�g���쐬�֐�
void mkdir(char dirname[]);
//TODO �����f�B���N�g���폜��A�f�B���N�g���쐬�֐� [�s�K�v]
//void delmkdirconfirm(char dirname[]);
//!? FEM�ɂ��CVCSEL���U���[�h�̓d���E���z�̎Z�o
void inputFEM();
//!? ���U���[�h�̑I��
void selectOsciMode(int Nvmode, double Vradius);	//TODO m, l, tau, beta, betaoverk, Rinf, eig, InfProf �������ɂ���		//char *directory
//!? -3dB�ш敝�v�Z�֐�
double sweep_g(int SingleSweep, double ginput, int Nvmode, int DMD, double r, double Vradius);

/* 1. �p�����[�^�̒�` */
// m: ���[�h����( TE&TM ~ 1, EH ~ n+1, HE ~ n-1), l: ���a�������[�h����, n: ���ʊp���[�h����
// N: ���a���Wr�̃R�A��������, Nbeta: �`���萔�̕�����, mater: �ޗ�ID, 
// lamda: �g��, A: �R�A���a, g: ���ܗ�����, n0: �R�A���S�̋��ܗ�, n1: �N���b�h�̋��ܗ�, dr: ���a���W���ݕ�
// k: �g��, delta: ����ܗ�, NA: �J����, aa: ���a���W�K�i���T�C�Y
// v: �K�i�����g��, w: �K�i���`���萔, D: �K�i���R�A�a
// tau: �Q���x, beta: �`���萔��, dbeta: �`���萔���ݕ�
// a: �W���s��S�̕��Ίp�v�f�i�[�z��, b: �W���s��S�̑Ίp�v�f�i�[�z��
// GI: �R�A�����ܗ�, q: �K�i���R�A�����ܗ�, q2: �R�A�����ܗ����U�p�����[�^
// R: �K�i���������d�ꐬ��, Rb: �t�ׂ���@�p����q�z��
// ML: LDU�����W���s���L�s�񕛑Ίp�v�f, MD: LDU�����W���s���D�s��Ίp�v�f
// dd, ds: �t�ׂ���@�ɂ�����������ʃp�����[�^, eig: �ŗL�l	
// modem: ���ʊp���[�h����, model: ���a���[�h����, modep: �僂�[�h����
// Nz: ���X�e�b�v��, Nzout: �t�@�C���o�͊�X�e�b�v��
// Li: i�Ԗڂ̔�����Ԃɂ����鑊�Βx���X�e�b�v���̍ő�l
// kim: i�Ԗڂ̔�����Ԃɂ�����e���[�h�̑��Βx���X�e�b�v��
// beta: �`���萔��, tau: �Q�x��, taumin: �ŏ��Q�x��, taumax: �ő�Q�x�� 
// zmax: �t�@�C�o��, dZ: ��ԍ��W���ݕ�, Tv: ���ԍ��ݕ�
// hmn: �d�͌����W��, Hmmmin: H�s��Ίp�v�f�ŏ��l, Hrowsum: 
// fmax: �ő�]�����g��, Nf: �]�����g���͈͂̕�����, df: �]�����g������\�i=fmax/Nf�j
// nmax: �v���t�@�C�����[�v���ɂ�����Lmax�̍ő�l�D
// nstd: ���Ύ��Ԃ̊�l�iy=0��tauimn����Ƃ���n=0�̕␳�l�j
// Dc: ���֒�, sigma2: microbending�̎��h�炬�̕��U
// wo: �ۓ��ݒ�p�p�����[�^�i0: w/o coupling, 1: w/ coupling_heterogeneity, 2: w/ coupling_microbending�j
// ftou: �t�@�C���o�͐ݒ�p�p�����[�^�i0: �����o��, 1: �S�o�́j
// matdis: �ޗ����U�l���p�����[�^�i0:����, 1:�l���j
// nP: �C���p���X�����i�[�z��p�m�ۗ̈�
// A00: ���̓C���p���X�U��, Aw00: ���̓C���p���X�U���̃X�y�N�g���������v
// �s�� A-(+)�i�e�ߓ_�ɂ����郂�[�h�����O��̃��[�h�p���[���z�j
// �s�� H�i�e�ߓ_�ɂ�����`�B�֐��j
// �z�� kim�ii�Ԗڂ̔�����Ԃɂ����郂�[�hm�̑��Βx���X�e�b�v���j

//! �ȉ����C���֐�   ////////////////////////////////////////////////////////////
// �P�F�̏ꍇ�̃��[�h���U�݂̂��l������EMBc���Z�o�������Ȃ�΁C
// OS_VCSEL_input_source_spectrum ��850nm�̂Ƃ��낾��1�ɑ���0�ɂ���
int main ( void ) {	
	
	FILE  *fp, *fptr, *fgvsBW, *frvsBW;
	char   directory[128];
	int    SingleSweep, rMeasure, DMD;
	double gmin, gmax, dg, ginput, bw;	 ginput = bw =0;
	double dtrsh; //�ǂݍ��ݗp
	double r = 0;
	double Vradius = 3.5;
	// double Vcore = 3.52;
	// double Vclad = 3.5;
	// double Vindexexp = 1.0e4;
	
	inputFEM();		//!? ��L��Vradius�Ȃǂ�ϐ��Ƃ��Ďg�p���Ă��悢���s�K�v�Ȃ̂ō���v����
	//sprintf(directory, "%s", inputFEM);

	//!? Nvmode, OffRes, OffRan �̓ǂݍ���
	int Nvmode;		double OffRes, OffRan;	
	if ((fp = fopen("[BW_Input].csv", "r")) != NULL) {
		char s1[128], s2[128];	int i = 0;
		for (int cnt=0 ; ; cnt++ ) {
			fscanf(fp, "%[^,], %[^,], ", s1, s2);
			if (strcmp(s2,"Nvmode")==0) { 
				fscanf(fp, "%d\n", &Nvmode);	i += 1; }
			else if (strcmp(s2, "DMD")==0) {
				fscanf(fp, "%d\n", &DMD); 		i += 1; }
			else if (strcmp(s2, "OffRes")==0) {
				fscanf(fp, "%lf\n", &OffRes);	i += 1; }
			else if (strcmp(s2, "OffRan")==0) {
				fscanf(fp, "%lf\n", &OffRan);	i += 1; }
			else { fscanf(fp, "%lfn", &dtrsh); }
			if (i == 4 || cnt > 200) { break; }		}}
	else { printf (" U cannot open the file !\n"); exit ( EXIT_FAILURE ); }
	fclose(fp);
	printf("Nvmode: %d\n", Nvmode);
	printf("OffRes: %.2lfum\n", OffRes);
	printf("OffRan: %.2lfum\n", OffRan);

	selectOsciMode(Nvmode, Vradius);		//!? m, l, tau, beta, betaoverk, Rinf, eig, InfProf �������ɂ���
	if ( (fptr = fopen("[BW_Input_index_exponent].csv","r")) != NULL ) {
		char ss1 [128], ss2 [128];
		 fscanf (fptr, "%[^,], %[^,], %d\n", ss1, ss2, &SingleSweep);
		 fscanf (fptr, "%[^,], %[^,], %lf\n", ss1, ss2, &gmin );
		 fscanf (fptr, "%[^,], %[^,], %lf\n", ss1, ss2, &gmax );
		 fscanf (fptr, "%[^,], %[^,], %lf\n", ss1, ss2, &dg );
		 //!? ���s�͍���폜����\��
		 fscanf(fptr, "%[^,], %[^,], %d\n", ss1, ss2, &rMeasure);	
		 fclose(fptr);}
	else { printf (" U cannot open the file !\n"); exit ( EXIT_FAILURE ); }

	//TODO �����ɃI�t�Z�b�g�ʂ��p�����[�^�Ƃ��� �֐�sweep_g �ɑ��
	if (SingleSweep == 0){
		if (DMD == 0){
			bw = sweep_g(SingleSweep, ginput, Nvmode, DMD, r, Vradius);		} //(�Ԃ�l��-3dB�ш敝)
		if (DMD == 1) {
			if ((frvsBW = fopen("BW_EMBc.csv", "r")) != NULL) {
				fprintf(frvsBW, "Offset r [um], -3dB bandwidth [GHz]\n");
				for (r = 0; r <= OffRan; r += OffRes) {
					bw = sweep_g(SingleSweep, ginput, Nvmode, DMD, r, Vradius);
					fprintf(frvsBW, "%.2lf, %lf\n", r, bw);		}}}}  //(�Ԃ�l��-3dB�ш敝)
	

	//! -3dB�ш敝�� g�l�ˑ���
	if (SingleSweep == 1){
		if ((fptr = fopen("[BW_Input_index_exponent].csv", "r")) != NULL) {
			if ((fgvsBW = fopen("BW_EMBc.csv", "w")) != NULL) {
				if (DMD == 0) {				
					fprintf(fgvsBW, "index exponent g, -3dB bandwidth [GHz]\n");
					for (int gcnt = 0; gcnt <= (gmax - gmin) / dg; gcnt++) {
						ginput = gmin + dg * gcnt;
						bw = sweep_g(SingleSweep, ginput, Nvmode, DMD, r, Vradius);  //(�Ԃ�l��-3dB�ш敝)
						fprintf(fgvsBW, "%.2lf, %lf\n", ginput, bw);	}}
				if (DMD == 1) {				
					fprintf(fgvsBW, "index exponent g, Offset r [um], -3dB bandwidth [GHz]\n");
					for (int gcnt = 0; gcnt <= (gmax - gmin) / dg; gcnt++) {
						ginput = gmin + dg * gcnt;
						for (r = 0; r <= OffRan; r += OffRes) {
							bw = sweep_g(SingleSweep, ginput, Nvmode, DMD, r, Vradius);  //(�Ԃ�l��-3dB�ш敝)
							fprintf(fgvsBW, "%.2lf, %.2lf, %lf\n", ginput, r, bw);	}}}
	fclose(fgvsBW);}}
	fclose(fptr);}		
}	
//! ���C���֐��I��   ////////////////////////////////////////////////////////////



//! �ȉ��T�u�֐�     ////////////////////////////////////////////////////////////
double sweep_g(int SingleSweep, double ginput, int Nvmode, int DMD, double r, double Vradius){
	//! //////////  �錾  //////////
	FILE  *fp,*fp2, *fp3, *fp4, *fq, *fr, *fr2, *fr3, *fr4, *fr5,*fr6, *fr7, *fr8, *fr9, *fr10, *fs, *fs2, *fs3, *fs4, *fdebug;
	int    i, j, l, m, n, y, nr, jmax, count = 0;
	int    mater, Nl, N, Nclad, Nbeta, NLP, NLP0, Ntotal, Nwkb, Ptotal, myu, nyu, nstd, nstdmin, nstdmax, mm;
	int    Nz, Nzout, Nf, Nfp, Ti, Li, Lmax, nmax, wo, gi, fout, matdis, scc, nP, Mn, launch, Nxy, Nom;
	double lamda, lamda0, lamdamin, lamdamax, lpmin, lpmax, dlp,  dl, k, omega,  A, AA, g, n0, n1, dr;
	double delta, NA, aa, v, w, D, w0, r0, dx, dy, xx, yy, rr, Rxy, Ein, Emev, Emod, Amev, Amod, r0dash;
	double tau, beta, dbeta, bb, eps1, eps2, sum, sumcore, sumclad, Rinf, dd, ds, eig;
	double deps2, Dc, sigma2, Db, dbmn, E_over,  hmn;
	double zmax, zout, dz, Tv, Hmmmin, Hrowsum, tauminstd, nctaumax, taumax, taumin = 0;
	double fmin, fmax, df, fpmin, fpmax, dfp, A00, Aw00, Hw, ReHw, ImHw, me, bw, tav, rms, Ptot, ap, spct, Cpsum;
	double *GI, *GC, *q, *q2, *R, *R2, *Rb, *a, *b, *ML, *MD, **Rlp, *Mtau, *Mbeta;
	double **H, *alpha, *P, *Pm, *Pg, *M, *Cp, *Pin, *Pout, *GIND, *OSin, *Wl, * WDin;
	int    *model, *modem,*modep, *pdeg, *kim, *Pnum;
	double V0, nv0, nv1, gsingle;	gsingle = 0;
	double dtrsh;
	double **OSmat;
	double ***Amin, ***Amplu;

	//! ////////////////////////////////////////////////////////////
	//! //////////            1-1. �����ݒ�                ///////////
	//! ////////////////////////////////////////////////////////////
	/*  ���̓t�@�C���̓ǂݍ��� */
	if ( (fp = fopen("[BW_Input].csv","r")) != NULL ) {
		char s1[128], s2[128];	
		 fscanf ( fp, "%[^,], %[^,], %d\n", s1, s2, &mater );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &gsingle);
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &A );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &zmax );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &lamda0 );
		 fscanf ( fp, "%[^,], %[^,], %d\n", s1, s2, &wo );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &deps2 );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &Dc );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &sigma2 );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &Db );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &dtrsh );
		 fscanf ( fp, "%[^,], %[^,], %d\n", s1, s2, &launch );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &r0dash);
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &dx );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &w0 );		
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &dtrsh);
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &dtrsh);
		 fscanf ( fp, "%[^,], %[^,], %d\n", s1, s2, &matdis );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &lamdamin );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &lamdamax );
		 fscanf ( fp, "%[^,], %[^,], %d\n", s1, s2, &Nl );
		 fscanf ( fp, "%[^,], %[^,], %d\n", s1, s2, &scc );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &lpmin );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &lpmax );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &dlp );
		 fscanf ( fp, "%[^,], %[^,], %d\n", s1, s2, &Ti );
		 fscanf ( fp, "%[^,], %[^,], %d\n", s1, s2, &gi );
		 fscanf ( fp, "%[^,], %[^,], %d\n", s1, s2, &fout );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &AA );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &dr );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &aa );
		 fscanf ( fp, "%[^,], %[^,], %d\n", s1, s2, &Nbeta );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &eps1 );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &eps2 );
		 fscanf ( fp, "%[^,], %[^,], %d\n", s1, s2, &jmax );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &dz );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &Tv );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &zout );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &fmin );
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &fmax );
		 fscanf ( fp, "%[^,], %[^,], %d\n", s1, s2, &Nf );
		 fscanf ( fp, "%[^,], %[^,], %d\n", s1, s2, &nP ); 
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &V0);
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &nv0);
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &nv1);	
		 fscanf ( fp, "%[^,], %[^,], %lf\n", s1, s2, &dtrsh);	}
	 else { printf (" U cannot open the file !\n"); exit ( EXIT_FAILURE ); }
	fclose(fp);

	//! g�l�̏㏑��
	 if (SingleSweep == 0){	g = gsingle;	}
	 if (SingleSweep == 1){	g = ginput;		}
	//! �I�t�Z�b�gr �̏㏑��
	 if (DMD == 0){	r0 = r0dash;}
	 if (DMD == 1){	r0 = r;		}

	 /* ���̓p�����[�^�̒P�ʕϊ� */
	 N = (int) ( 1000.0*A / dr ); Nclad = (int)( 1000.0*( AA - A ) / dr ); // �ʓ����a�����X�e�b�v���̊��Z�i�K�����̈ʒu�Œ�`�j�D
	 if ( Nl%2 != 0 ) { Nl = Nl+1;} dl = ( lamdamax - lamdamin ) / (double)Nl;  // �ޗ����U�]���p�g���X�e�b�v	
	 fpmin = C / lpmax, fpmax = C / lpmin, Nfp = int (( lpmax - lpmin ) / dlp), dfp = (fpmax-fpmin) / (double) Nfp; // ���g���̈�̌����X�y�N�g���̒�`�D
	 lamda0 = lamda0*1.0e-9; lamdamin = lamdamin*1.0e-9; lamdamax = lamdamax*1.0e-9; dl = dl*1.0e-9; // �P�ʕϊ��im�j
	 lpmin = lpmin*1.0e-9, lpmax = lpmax*1.0e-9, dlp = dlp*1.0e-9; // �P�ʕϊ��im�j
	 dfp = dfp*1.0e-3; fpmin = fpmin*1.0e-3; fpmax = fpmax*1.0e-3; // �P�ʕϊ��iTHz�j
	 A = A*1.0e-6; dr = dr*1.0e-9; AA = AA *1.0e-6; // �P�ʕϊ��im�j
	 w0 = w0*1.0e-6; r0 = r0*1.0e-6, dx = dy = dx*1.0e-9; // �P�ʕϊ��im�j
	 Nxy = (int) ( 5.0*w0 / dx );//�Ƃ肠����4w0�͈̔�
	 Dc = Dc*1.0e-9, Db = Db*1.0e-3, deps2 = deps2 *(Epsilon0)*(Epsilon0); // �P�ʕϊ��im�j�C��U�d����U�d���ɕϊ�
	 fmin = fmin*1.0e-3; fmax = fmax*1.0e-3; df = ( fmax-fmin ) / (double) Nf;	// �P�ʕϊ��iTHz(1/ps)�j
	 Nz = (int) ( zmax / dz ); Nzout = (int) ( zout / dz ); // �t�@�C�o����������������уt�@�C���o�͊Ԋu�̊��Z�D
	 D = A / aa; // �K�i���R�A�a�̊��Z�D
	 lamda = lamda0; Aw00 = 0.0; spct = Cpsum = 0.0; // �ϐ�������
	 xx = yy = rr = Rxy = 0.0;
	 if ( matdis == 0 ) { Nl = 0; }
	
	 /* ���̓f�[�^�̊i�[ */
	 // ���ˌ����Ԕg�`
	 Pin = drealvector ( 0, Ti ); init_realvector ( Pin, 0, Ti );
	 if ( (fr5 = fopen ("Input_pulse_waveform.csv", "r") ) != NULL ) {
		 for ( i = 0; i < Ti; i++ ) { fscanf ( fr5, "%lf\n", &Pin[i] ); }}
	 else { printf (" U cannot open the file !\n"); exit ( EXIT_FAILURE );	}
	 fclose (fr5);
	 // ���ˌ��X�y�N�g��
	 
	 OSin = drealvector ( 0, Nfp );	init_realvector ( OSin, 0, Nfp );
	 Wl = drealvector ( 0, Nl+1 );	init_realvector ( Wl, 0, Nl+1 );
	 OSmat = dmatrix(0, Nfp, 0, Nvmode);	init_dmatrix(OSmat, 0, Nfp, 0, Nvmode);
	 //���˃X�y�N�g���̓ǂݍ���
	 if (launch == 0 || launch == 1) {
		 Nvmode = 1;
		 if ((fr8 = fopen("[OS_Input_source_spectrum].csv", "r")) != NULL) {
			 for (i = 0; i <= Nfp; i++) { j=0;
				 //for (j = 0; j < Nvmode; j++) {
					 fscanf(fr8, "%lf\n", &OSmat[i][j]);	}}//}
		 else { printf(" U cannot open the file !\n"); exit(EXIT_FAILURE); }
		 fclose(fr8);
	 }
	 if (launch == 2) {
		 if ( (fr8 = fopen ("[OS_VCSEL_input_source_spectrum].csv", "r") ) != NULL ) {
			 for ( i = 0; i <= Nfp; i++ ) { 
				 for (j = 0; j < Nvmode; j++) {
					 fscanf ( fr8, "%lf\n", &OSmat[i][j] ); }}}//, %lf	&Wl[i],
		 else { printf (" U cannot open the file !\n"); exit ( EXIT_FAILURE );	}
		 fclose(fr8);
	 }
	 // ���ˌ��X�y�N�g���ɂ����锭�U���[�h�̔��ʔg�������
	 if ( (fr10 = fopen ("[OS_Wavelength_discrimination].csv", "r") ) != NULL ) {	//TODO WDin[ ]�͋��E�g��������
		 char s1[128], s2[128];
		 fscanf(fr10, "%[^,], %[^,], %d\n", &s1, &s2, &Nom);	printf("%d\n", Nom);
		 WDin = drealvector(0, Nom);	init_realvector(WDin, 0, Nom);
		 for ( i = 0; i < Nom; i++ ) { fscanf ( fr10, "%lf\n", &WDin[i] );	 printf("%lf\n", WDin[i]);}}//, %lf	&Wl[i],
	 else { printf (" U cannot open the file !\n"); exit ( EXIT_FAILURE );	}
	 fclose(fr10);
	 // ���ܗ����z�i@589nm�j&�h�[�p���g�Z�x���z���z
	 if ( gi == 1) {	 
		 GIND = drealvector ( 0, N ); init_realvector ( GIND, 0, N );	 
		 if ( (fr7 = fopen ("GI_profile_NaD.csv", "r") ) != NULL ) {
			 for ( i = 0; i <= N; i++ ) { fscanf ( fr7, "%lf\n", &GIND[i] ); }}
		 else { printf (" U cannot open the file !\n"); exit ( EXIT_FAILURE );	}
		 GC = drealvector ( 0, N ); init_realvector ( GC, 0, N );	 
		 for ( i = 0; i <= N; i++ ) { GC[i] = ( GIND[i] - GIND[N] ) / dndc ( 589.0, mater); }	
		 fclose (fr7);}

	 //! ���Ԃ�����Έȉ����s������
	 //  ���Ԃ�����΁A�o�̓f�B���N�g���̐ݒ���s������
	 //  (�Ȃ��Ȃ��f�o�b�O�ł��Ȃ��ꍇ�͉��L���쐬����K�v����)
	 /*	 fopen�֐���write�^�C�v�ŊJ���t�@�C���́A���L�̏o�̓t�@�C���̐ݒ�őS�Ăł��邽�߁A
		 ���L���f�B���N�g���Ɏ��߂�
		 g=%.2lf���܂߂Ċ֐��𕪂��A�Ԃ�l -3dB�ш敝*�ł�肽��	
		 �ҏW����Ȃ�A
		 �@�iSingleSweep==1�j�f�B���N�g���Eg�l��2��BW_setting�Ȃǂ̏o�̓t�@�C�����ɂ���iex. BW_setting\�j
		 �A�iSingleSweep==0�j��コ�񂩂���������Ԃ̂܂܏o�̓t�@�C���̐ݒ���s��*/

	 //! �o�̓t�@�C���̐ݒ�
	 // �]�������ꗗ
	 //! sprintf(frpath, "%s/profile.csv", directory);
	 if ( ( fp2 = fopen("BW_setting.csv","w")) != NULL ) {}
	 else { printf (" U cannot open the file !\n"); exit ( EXIT_FAILURE ); }
	 // ���ܗ��v���t�@�C��
	 if ( ( fp3 = fopen ("BW_profile.csv", "w")) != NULL ) {
		 if (gi == 0) { fprintf (fp3, "r [um],n,q,q2\n"); } if (gi == 1) { fprintf (fp3, "r [um],c [wt pct],n,q,q2\n"); }}	
	 else { printf (" U cannot open the file !\n"); exit ( EXIT_FAILURE ); } 
	 // �L���v�f�@�v�Z����
	 if ( ( fp4 = fopen ("FEM_result.csv", "w")) != NULL ) {}	  			
	 else { printf (" U cannot open the file !\n"); exit ( EXIT_FAILURE ); }
	 // �d�͌����W���v�Z����
	 if ( ( fq = fopen ( "Hmn_result.csv", "w" )) != NULL ) {
		 if ( wo == 1 ) {			
			 fprintf ( fq, "��,m (mode ��),l (mode ��),p (mode ��),�� (mode ��),��,m (mode ��),l (mode ��),p (mode ��),�� (mode ��),��p,|����|,E-field overlap,h�ʃ�\n");}
		 if ( wo == 2 ) {
			 fprintf ( fq, "��,m (mode ��),l (mode ��),p (mode ��),�� (mode ��),��,m (mode ��),l (mode ��),p (mode ��),�� (mode ��),��p,|����|,h�ʃ�\n");}}				
	 else { printf ( "U cannot open the file!\n" ); exit ( EXIT_FAILURE ); }
	 // �s��H
	 if ( ( fr = fopen ( "CPE_Hmatrix.csv", "w" ) ) != NULL ) {}
	 else { printf ( "The file cannot be opened !\n" ); exit ( EXIT_FAILURE ); }			
	 // �C���p���X�����g������ P
	 if ( ( fr2 = fopen ( "CPE_Impulse-responce.csv", "w" ) ) != NULL ) {}			
	 else { printf ("The file cannot be opened !\n"); exit ( EXIT_FAILURE ); }
	 // �o�˔g�`�g������ Pout
	 if ( ( fr6 = fopen ( "CPE_Output-pulse-waveform.csv", "w" ) ) != NULL ) {}			
	 else { printf ("The file cannot be opened !\n"); exit ( EXIT_FAILURE ); }
	 // ���[�h�p���[���z Pm		
	 if ( ( fr3 = fopen ( "CPE_Mode-power-distribution.csv", "w" ) ) != NULL ) {}			
	 else { printf ("The file cannot be opened !\n"); exit ( EXIT_FAILURE ); }
	 // ���[�h�Q�p���[���z Pg		
	 if ( ( fr9 = fopen ( "CPE_Group-power-distribution.csv", "w" ) ) != NULL ) {}			
	 else { printf ("The file cannot be opened !\n"); exit ( EXIT_FAILURE ); }	
	 // ���g������ M
	 if ( ( fr4 = fopen ( "CPE_Frequency-respnse.csv", "w" ) ) != NULL ) {
		 fprintf ( fr4, "," ); for ( j = 0; j <= Nf; j++ ) { fprintf ( fr4, "%f,", (fmin*1.0e3)+(double)j*(df*1.0e3) ); } fprintf ( fr4, "\n" ); }
	 else { printf ("The file cannot be opened !\n"); exit ( EXIT_FAILURE ); }
	 // ��v�Ȍ���
	 if ( ( fs = fopen ( "BW_result.csv", "w" )) != NULL ) {
		 if ( matdis == 0) { fprintf ( fs, "g value,Nlp,Nlp0,tstart [ns],length[m],pulse broadening [ps],-3dB bandwidth [GHz],tav[ps],Ptot,Aw00,spct\n"); }
		 if ( matdis == 1) { fprintf ( fs, "g value,length[m],pulse broadening [ps],-3dB bandwidth [GHz],tav[ps],Ptot,Aw00,spct\n"); }}
//		 else { fprintf ( fs, "g value,length[m],-3dB bandwidth [GHz],rms width [ps],tav[ps],Ptot,Aw00,spct\n"); }}
	 else { printf ( "U cannot open the file!\n" ); exit ( EXIT_FAILURE ); }
	 // �C���p���X����
	 if ( ( fs3 = fopen ( "BW_Impulse-responce.csv", "w" )) != NULL ) {
		 for ( n = 0; n <= nP; n++ ) { fprintf ( fs3, "%f,", double(n)*Tv ); } fprintf ( fs3, "\n" ); }
	 else { printf ( "U cannot open the file!\n" ); exit ( EXIT_FAILURE ); }
	 // �o�˔g�`�g������ Pout
	 if ( ( fs4 = fopen ( "BW_Output-pulse-waveform.csv", "w" ) ) != NULL ) {}			
	 else { printf ("The file cannot be opened !\n"); exit ( EXIT_FAILURE ); }
	 // �����X�y�N�g��
	 if ( ( fs2 = fopen ( "BW_source-spectrum.csv", "w" )) != NULL ) {
			 for ( i = 0; i <= Nfp; i++ ) { fprintf ( fs2, "%f,%f,%f\n", fpmin+double(i)*dfp,1.0e-3*C / ( fpmin+double(i)*dfp ), OSin[Nfp-i] ); }}
	 else { printf ( "U cannot open the file!\n" ); exit ( EXIT_FAILURE ); }
	 if ( ( fdebug = fopen ( "[Debug].csv", "w" )) != NULL ) {
		 fprintf ( fs2, "Debug start\n");}
	 else { printf ( "U cannot open the file!\n" ); exit ( EXIT_FAILURE ); }

	 printf ( "Total spatial steps Nf: %d\n", Nf );

	 //! ////////////////////////////////////////////////////////////
	 //! //////////           1-2. Input FEM              ///////////
	 //! ////////////////////////////////////////////////////////////

	 if (launch != 2) {
		 goto SkipInputFEM;	 }

 SkipInputFEM:

		//! //////////////////////////////////////////////////////////////////
		//! //////////      ( 2. ~ 5. �̃��[�v�ɂ�����ݒ�̏o�� )        //////////
		//! /////////////////////////////////////////////////////////////////
		/* �ϐ��̏����� */
		Aw00 = tauminstd = 0.0; nmax = nstd = nstdmin = nstdmax = Mn = 0;
		/* �z��̋L���̈�m�ۂ���я����� */
		if ( matdis == 1 ) { P = drealvector ( 0, nP ), init_realvector ( P, 0, nP ); }
		Pnum= dintvector ( 0, Nl ); init_intvector ( Pnum, 0, Nl );
		M = drealvector ( 0, Nf ); init_realvector ( M, 0, Nf );
		int contct=0;
		for ( y = 0; y <= Nl; y++ ) {					   
			/* �g���w�� */
			if ( matdis == 1 ) { lamda = lamdamax -  y*dl; }
			/* �e��p�����^�v�Z */
			if ( gi == 0 ) { n1 = ncore ( lamda*1.0e9, mater ); n0 = nclad ( lamda*1.0e9, mater ); }
			if ( gi == 1 ) { n1 = dndc ( lamda*1.0e9, mater )*GC[0] + nclad ( lamda*1.0e9, mater ); n0 = nclad ( lamda*1.0e9, mater ); }	
			delta = (n1*n1-n0*n0) / (2.0*n1*n1), NA = sqrt (n1*n1-n0*n0);
			k = 2.0*PI / lamda; omega = 2.0*PI*C / lamda; v = k*aa*n1*sqrt ( 2.0*delta );
			Nwkb = int( (1.0/4.0)*( g / (g+2.0) )*(k*k)*(n1*n1)*delta*(A*A) );
			/* �ϐ��̏����� */
			bb = -1; dbeta=  k*( n1 - n0 ) / (double) Nbeta;
			NLP= 0;  NLP0 = Ntotal = Ptotal = 0;
			/* �z��̋L���̈�m�ۂ���я����� */
			GI = drealvector ( 0, N ); init_realvector ( GI, 0, N );
			q = drealvector ( 0, N ); init_realvector ( q, 0, N );
			q2 = drealvector ( 0, N ); init_realvector ( q2, 0, N );
			R= drealvector ( 0, N ); init_realvector ( R, 0, N );
			R2 = drealvector ( 0, N ); init_realvector ( R2, 0, N );
			Rb = drealvector ( 0, N ); init_realvector ( Rb, 0, N );
			a = drealvector ( 0, N ); init_realvector ( a, 0, N );
			b = drealvector ( 0, N ); init_realvector ( b, 0, N );
			ML = drealvector ( 0, N ); init_realvector ( ML, 0, N );
			MD = drealvector ( 0, N ); init_realvector ( MD, 0, N );
			modem = dintvector ( 0, 2*Nwkb ); init_intvector ( modem, 0, 2*Nwkb );
			model = dintvector ( 0, 2*Nwkb ); init_intvector ( model, 0, 2*Nwkb );
			modep = dintvector ( 0, 2*Nwkb ); init_intvector ( model, 0, 2*Nwkb );
			pdeg = dintvector ( 0, Nwkb ); init_intvector ( pdeg, 0, Nwkb );
			Mbeta = drealvector ( 0, 2*Nwkb ); init_realvector ( Mbeta, 0, 2*Nwkb );
			Mtau = drealvector ( 0, 2*Nwkb ); init_realvector ( Mtau, 0, 2*Nwkb );
			Rlp = dmatrix ( 0, 2*Nwkb, 0, N ); init_dmatrix ( Rlp, 0, 2*Nwkb, 0, N );
			/* ���ܗ����z�C����ܗ����z����єg���������z */
			if (gi == 0) {		
				for ( j = 0; j <= N; j++ ) { GI[j] = n1*sqrt ( 1.0-2.0*delta*pow (((double) j / (double) N), g) ); } 
				for ( j = 0; j <= N; j++ ) { q[j] = (GI[j]*GI[j] - n0*n0) / (n1*n1 - n0*n0); }		
				for ( j = 0; j <= N; j++ ) { q2[j] = GI[j]*GI[j] - (lamda*GI[j]*dndl (lamda*1.0e9, GI[j], mater)) / (1 - (lamda/GI[j])*dndl (lamda*1.0e9, GI[j], mater)); }				
				if ( y == Nl/2 || fout == 1 ) { for ( j = 0; j <= N; j++ ) { fprintf ( fp3, "%f, %f, %f, %f\n", A*1.0e6* (double) j / (double) N, GI[j], q[j], q2[j] ); }}}
			if (gi == 1) {
				for ( j = 0; j <= N; j++ ) { GI[j] = dndc ( lamda*1.0e9, mater )*GC[j] + nclad ( lamda*1.0e9, mater ); } // �]���g���ɂ�������ܗ����z�Ɋ��Z
				for ( j = 0; j <= N; j++ ) { q[j] = (GI[j]*GI[j] - GI[N]*GI[N]) / (GI[0]*GI[0] - GI[N]*GI[N]); }
				for ( j = 0; j <= N; j++ ) { q2[j] = GI[j]*GI[j] - (lamda*GI[j]*dndl (lamda*1.0e9, GI[j], mater)) / (1 - (lamda/GI[j])*dndl (lamda*1.0e9, GI[j], mater)); }
				if ( y == Nl/2 || fout == 1 ) { for ( j = 0; j <= N; j++ ) { fprintf ( fp3, "%f, %f, %f, %f, %f\n", A*1.0e6* (double) j / (double) N, GC[j], GI[j], q[j], q2[j] ); }}}

			//! �]�������̏o�� 
			if (  y == 0 ) {
				fprintf ( fp2, "Material, mater, %d\n", mater );
				fprintf ( fp2, "Index exponent, g, %.2lf\n", gsingle );
				fprintf ( fp2, "Central wavelength, ��0, %f nm\n", lamda0*1e9 );
				fprintf ( fp2, "Step size of wavelength, dl, %f nm\n", dl*1e9 );
				fprintf ( fp2, "Partition number of evaluated wavelength, Nl, %d\n", Nl );
				fprintf ( fp2, "Core radius, A, %f um\n", A*1e6 );
				fprintf ( fp2, "Analysis region in radial axis, AA,%f um\n", AA*1e6 );
				fprintf ( fp2, "Refractive index at the core center, n1,%f\n", n1 );
				fprintf ( fp2, "Refractive index in the cladding, n0,%f\n", n0 );
				fprintf ( fp2, "Relative refractive index, ��,%f\n", delta );
				fprintf ( fp2, "Numerical aperture, NA, %f\n", NA );
				fprintf ( fp2, "Step size of the elements, dr, %f, nm\n", dr*1e9 );
				fprintf ( fp2, "Step size of propagation constant, d��,%f\n", dbeta );
				fprintf ( fp2, "Partition number of propagation constant, N��,%d\n", Nbeta );
				fprintf ( fp2, "Partition number of fiber core radius, N,%d\n", N );
				fprintf ( fp2, "Partition number of fiber cladding, Nclad,%d\n", Nclad );
				fprintf ( fp2, "Maximum allowable error for convergence solution vector, eps1,%f\n", eps1 );
				fprintf ( fp2, "Maximum allowable error of zero eigen value, eps2,%f\n", eps2 );
				fprintf ( fp2, "Maximum number of iterations in inverse power method, jmax,%d\n", jmax );
				fprintf ( fp2, "Correlation length, Dc, %e, nm\n", Dc*1e9 );
				fprintf ( fp2, "Mean square of dielectric constant fluctuation, <d��2>, %e\n", deps2 );
				fprintf ( fp2, "Total fiber length,zmax,%e,m\n", zmax );
				fprintf ( fp2, "Step size of fiber length,dz,%e,m\n", dz );
				fprintf ( fp2, "Step size of time,Tv,%e,ps\n", Tv );
				fprintf ( fp2, "Total spatial steps,Nz,%d\n", Nz );
				fprintf ( fp2, "File output step interval,Nzout,%d\n", Nzout );
				fprintf ( fp2, "Minimum evaluated frequency,fminx,%e,GHz\n", fmin*1.0e3 );
				fprintf ( fp2, "Maximum evaluated frequency,fmax,%e,GHz\n", fmax*1.0e3 );
				fprintf ( fp2, "Step size of evaluated frequency,df,%e,GHz\n", df*1.0e3 );	
				printf ( "Central wavelength ��0: %f nm\n", lamda0*1e9 );
				printf ( "Step size of wavelength dl: %f nm\n", dl*1e9 );
				printf ( "Partition number of evaluated wavelength Nl: %d\n", Nl );
				printf ( "Minimum evaluated spectral frequency, fpmin, %f THz\n", fpmin );
				printf ( "Maximum evaluated spectral frequency, fpmax, %f THz\n", fpmax );
				printf ( "Core radius A: %f um\n", A*1.0e6 );
				printf ( "Analysis region AA: %f um\n", AA*1e6 );
				printf ( " Refractive index at the core center n1: %f\n", n1 );
				printf ( "Refractive index in the cladding n0: %f\n", n0 );
				printf ( "Relative refractive index ��: %f\n", delta );
				printf ( "Numerical aperture NA: %f\n", NA );
				printf ( "Step size of the elements dr: %f nm\n", dr*1.0e9 );
				printf ( "Partition number of fiber core radius N: %d\n", N );
				printf ( "Partition number of fiber cladding Nclad: %d\n", Nclad );
				printf ( "Step size of propagation constants d��: %f\n", dbeta );
				printf ( "Partition number of propagation constants N��: %d\n", Nbeta );
				printf ( "Maximum allowable error for convergence solution vector eps1: %f\n", eps1 );
				printf ( "Maximum allowable error of zero eigen value eps2: %f\n", eps2 );
				printf ( "Maximum number of iterations in inverse power method jmax: %d\n", jmax );
				printf ( "Correlation length Dc: %e m\n", Dc );
				printf ( "Mean square of dielectric constant fluctuation <d��2>: %f\n", deps2 );
				printf ( "Total fiber length zmax: %e m\n", zmax );
				printf ( "Step size of fiber length dz: %e m\n", dz );
				printf ( "Step size of time Tv: %e ps\n", Tv );
				printf ( "Total spatial steps Nz: %d\n", Nz );
				printf ( "File output step interval Nzout: %d\n", Nzout );
				printf ( "Minimum evaluated frequency fminx: %e GHz\n", fmin*1.0e3 );
				printf ( "Maximum evaluated frequency fmax: %e GHz\n", fmax*1.0e3 );
				printf ( "Step size of evaluated frequency df: %e GHz\n\n", df*1.0e3 ); }				
			printf ( "Wavelength: %f nm (%d/%d)\n", lamda*1.0e9,y,Nl );
			if ( wo ==0 ) { printf ( "without mode coupling\n"); }
			if ( wo ==1 ) { printf ( "with microscopic heterogeneities\n"); }
			if ( wo ==2 ) { printf ( "with microbending\n"); }
			fprintf ( fp2, "lamda=%f nm\n", lamda*1.0e9 );
			fprintf ( fp4, "lamda=%f nm\n", lamda*1.0e9 );
			if (  y == Nl/2 || fout == 1 ) {
				fprintf ( fp4, "m,l,p,tau,beta,ne,confinement,R_infinite,eig,");
				for ( j = 0; j <= N; j++ )  { fprintf ( fp4,"%f,", (double) j*dr*1.0e6 ); } fprintf ( fp4,"\n" ); }
			fprintf ( fq, "lamda=%f nm\n", lamda*1.0e9 );
			fprintf ( fr2, "lamda=%f nm\n", lamda*1.0e9 );
			fprintf ( fr3, "lamda=%f nm\n", lamda*1.0e9 );


			//! /////////////////////////////////////////////////////////////////
			//! //////////           2. ���[�h��� (FEM.cpp)             //////////
			//! ////////////////////////////////////////////////////////////////
			//
			for ( m = 0 ; ; m++ ) {
				l = 1; beta = k*n1;
				/****************************************************************/
				for ( ; ; )  {
					/* �W���s��v�Z����щ����R���X�L�[���� */		
					w = aa*sqrt ( (beta*beta) - (k*k)*(n0*n0) );
					S_matrix ( a, b, q, m, N, v, w, D );
					mcholesky ( a, b, ML, MD, m, N );
					/* �����x�N�g�� R0 �̕t�^ */
					R0 ( MD, R, m, N );
					/*  �A���ꎟ������ SR=(LDL)R=bR �̔����]�� */
					for ( j = 0 ; ; j++ ) { 
						for ( i = 0; i <= N; i++ )  { Rb[i] = R[i]; }
						mcholesky_sol ( ML, MD, R, m, N );
						R_norm ( R, N ); 
						/* �������� */
						dd = 0, ds = 0;
						for ( i = 0; i <= N; i++ ) {
							dd = dd + (Rb[i] - R[i])*(Rb[i] - R[i]);
							ds = ds + (Rb[i] + R[i])*(Rb[i] + R[i]); }
						if ( dd < eps1 || ds < eps1 ) break;
						if ( j >= jmax ) goto next;
						// �@ R��Rb�̐�����dd��0�ɑQ�߂���Ύ����ibreak�j�D
						// �A R��Rb�̐����ads��0�ɑQ�߂���Β��~�ibreak�j�D
						// �B �����񐔂�����l jmax �𒴂��������ύX���čČv�Z�igo to next�j�D
					}
					/* �ŗL�l�̌v�Z */
					eig = Eigen ( R, a, b, N, m );
					/* �ŗL�l�̑Ó����]�� */
					// �u�����ŗL�l eig ���O��l bb �Ɠ��l�v�ł���΃���ς��ď��߂���Čv�Z
					if ( eig == bb ) {
						dbeta = k*(n1 - n0) / (double) Nbeta;
						beta = beta - 1.0*dbeta;
						continue; }
				 
					//  �@�u0< �����ŗL�l eig < eps2�v�ł���Η�ŗL�l�Ƃ��č̗p
					if ( 0.0 < eig && eig < eps2 )  {
						/* �������d�ꐬ��R�̋K�i�� �i�p���[��1W�Ƃ���j*/			  
						sum = sumcore = sumclad = 0.0;
						for ( j = 0; j < N; j++ ) { sumcore = sumcore + R[j]*R[j]*(j*dr)*dr; }
						for ( j = 0; j < Nclad; j++ ) { sumclad = sumclad + R[N]*( bessk (m, w*(j*dr+A)) / bessk (m, w*A) )*R[N]*( bessk (m, w*(j*dr+A)) / bessk (m, w*A) )*(j*dr+A)*dr; }
						for ( j = 0; j <= N; j++ ) { R2[j] = R[j] * sqrt ( (omega*Mu0) / (PI*beta*(sumcore + sumclad)) ); }
						tau = (1.0 / (C*1.0e-12))*(k / beta) * dbdk_bunshi ( R, q2, D, w, m, N ) / dbdk_bunbo ( R, D, w, m, N ); // = (1/c)*(d��/dk) [ps/m]
						/* �v�Z���ʂ̏o�� */
						Rinf = R[N]*( bessk (m, w*(Nclad*dr+A)) / bessk (m, w*A) );	
						if ( y == Nl/2 || fout == 1 ) { fprintf ( fp4, "%d, %d, %d, %f, %f, %f, %f, %f, %f,",m,l,2*l+m-1,tau,beta,beta/k,sumcore/(sumcore+sumclad),Rinf,eig ); }
						modem[NLP] = m, model[NLP] = l, Mtau[NLP] = tau, Mbeta[NLP] = beta;
						modep[NLP] = 2*l + m - 1; if ( Ptotal < modep[NLP] ) { Ptotal = modep[NLP]; } // �僂�[�h����
						for ( j = 0; j <= N; j++ )  { if ( fout == 1 || y == Nl/2 ) { fprintf ( fp4,"%f,", R2[j] ); } Rlp[NLP][j] = R2[j]; }
						if ( y == Nl/2 || fout == 1 ) { fprintf ( fp4,"\n" ); }
						//printf ( "%d, %d, %f, %f, %f, %f, %f\n", m, l, tau, beta, beta/k, eig, R[N] );

						dbeta = k*(n1 - n0) / (double) Nbeta;
						bb = eig;
						l = l + 1;
						NLP = NLP + 1; if ( m == 0 ) { NLP0 = NLP0 +1;}
						count = 0; }

					//  �A�u0 < �����ŗL�l eig�v���u-1 < �O��l bb < 0�v�ł����
					if ( eig > 0.0 && bb < 0.0 && (fabs(bb) < 1.0) ) {
						beta = beta + dbeta;
						dbeta = dbeta / 2.0;
						count = count + 1; }

					// �B ���̑�
					else { bb = eig; count = 0; }

					if ( count > 1000 ) {	
						dbeta = k*(n1 - n0) / (double) Nbeta;
						beta = beta - dbeta; }
next:
					beta = beta - dbeta;
					if ( beta < k*n0 )	break; 
				}
				/****************************************************************/
				if ( l == 1 ) { Ntotal = ( NLP0 )*2 + ( NLP - NLP0 )*4; break; } // �����ō������ɓ��B
			}

			free_drealvector ( GI, 0 ); free_drealvector ( q, 0 ); free_drealvector ( q2, 0 );
			free_drealvector ( R, 0 ); free_drealvector ( R2, 0 ); free_drealvector ( Rb, 0 );
			free_drealvector ( a, 0 ); free_drealvector ( b, 0 );
			free_drealvector ( ML, 0 );	free_drealvector ( MD, 0 );

			printf ( "Total numbers of LPml modes (WKB) Nwkb: %d\n", Nwkb);
			printf ( "Total numbers of LPml modes NLP: %d\n", NLP);
			printf ( "Total numbers of LP0l modes NLP0: %d\n", NLP0);
			printf ( "Total numbers of all modes Ntotal: %d\n", Ntotal); 
			printf ( "Total numbers of all mode groups Ptotal: %d\n", Ptotal); 


			//! ///////////////////////////////////////////////////////////////
			//! //////////           3. LP���[�h�����̐���              ///////////
			//! ///////////////////////////////////////////////////////////////
			/* �Q�x���͈� */
			// �e���[�v�̒P�ʒ��Q�x���͈́itaumin < tau < taumax�j
			taumax = taumin = Mtau[0];
			for ( myu = 0; myu < NLP; myu++ ) {
				if ( taumax < Mtau[myu] ) { taumax = Mtau[myu]; }
				if ( taumin > Mtau[myu] ) { taumin = Mtau[myu]; }}
			// �e���[�v�̍ō������[�h�Q�P�ʒ��Q�x��
			nctaumax = taumin;
			for ( myu = 0; myu < NLP; myu++ ) {							
				if ( modep [myu] == Ptotal ) {	if ( nctaumax < Mtau[myu] ) { nctaumax = Mtau[myu]; }}}
			// �e���[�v�̊���ԕ␳�v�f���instd�j
			if ( y == 0 ) { tauminstd = taumin; }
			nstd = (int) ( (taumin - tauminstd)*zmax / Tv );
			printf ("nstd:%d,nstdmin:%d\n", nstd,nstdmin);
			fprintf ( fp2, "Minimum group delay per unit length,taumin,%e,ps/m\n", taumin );
			fprintf ( fp2, "Maximum group delay per unit length,taumax,%e,ps/m\n", taumax );
			// �e���[�v�̍ő�S�x�����iLmax�j
			Lmax = (int) ( (( taumax - taumin )*( (double)Nz*dz ) / Tv) + 0.5 );
			fprintf ( fp2, "Total time steps,Lmax,%d\n", Lmax );
			if ( Lmax > ( 1.0e9 / (8.0*NLP) ) ) { printf ( "Memory over!\n" ); exit ( EXIT_FAILURE ); }
			printf ( "Required memory for CPE analysis: %fGB\n", (double) (Lmax*NLP*8) /1.0e9 );
			// �v���t�@�C�����[�v�̍ő�S�x�����inmax�j
			if ( Lmax > nmax ) { nmax = Lmax; }
			fprintf ( fp2, "Maximum time step,nmax,%d\n", nmax );
			// �C���p���X�����i�[�z�񐔁iPnum�j
			if ( nstd == 0 ) { Pnum [y] = nmax + ( nstdmax - nstdmin ); }
			if ( nstd > 0 ) {
				if ( nstd < nstdmax ) { Pnum [y] = nmax + ( nstdmax - nstdmin ); }
				else { Pnum [y] = nmax + ( nstdmax - nstdmin ) + ( nstd - nstdmax ); }}
			if ( nstd < 0 ) {
				if ( nstd < nstdmin ) { Pnum [y] = nmax + ( nstdmax - nstdmin ) + ( nstdmin - nstd ); }
				else { Pnum [y] = nmax + ( nstdmax - nstdmin ); }}
			fprintf ( fp2, "Net total time steps,Pnum,%d\n", Pnum [y] );

			/* �z��̋L���̈�m�ۂ���я����� */	
			H = dmatrix ( 0, NLP-1, 0, NLP-1 );			init_dmatrix ( H, 0, NLP-1, 0, NLP-1 );
			Amin = d3tensor ( 0, Nvmode, 0, NLP-1, 0, Lmax );
			init_d3tensor ( Amin, 0, Nvmode, 0, NLP-1, 0, Lmax );
			alpha = drealvector ( 0, NLP-1 );			init_realvector ( alpha, 0, NLP-1 );
			kim = dintvector ( 0, NLP-1 );				init_intvector ( kim, 0, NLP-1 );
			Pm = drealvector ( 0, NLP-1 );				init_realvector ( Pm, 0, NLP-1 );
			Pg = drealvector ( 0, Ptotal-1 );			init_realvector ( Pg, 0, Ptotal-1 );
			Amplu = d3tensor(0, Nvmode, 0, NLP - 1, 0, Lmax);	
			init_d3tensor(Amplu, 0, Nvmode, 0, NLP - 1, 0, Lmax);
			
			if ( matdis == 0 ) { P = drealvector ( 0, Lmax ); init_realvector ( P, 0, Lmax ); }
			A00 = 0.0;		

			//! /////////////////////////////////////////////////////////////////////////
			//! //////////      4. LP���[�h �� ��LP���[�h �� �Ԃ̓d�͌����W���v�Z        //////////
			//! ////////////////////////////////////////////////////////////////////////
			//!
			//! H�s��̎Z�o
			if ( wo == 0 ) { for ( myu = 0; myu < NLP; myu++ ) { H[myu][myu] = 1.0; } }
			else {
			// �~�N���s�ψ�\��
			if ( wo == 1) { dbmn = 0.0; E_over = 0.0; hmn =0.0;
#pragma omp parallel			
			{
#pragma omp for
			for ( myu = 0; myu < NLP; myu++ ) {							
				for ( nyu = 0; nyu < NLP; nyu++ ) {
					dbmn = Mbeta [myu] - Mbeta [nyu];						
					for ( n = 1; n <= N; n++ ) {												
						E_over = E_over + (2.0*PI)*( Rlp[myu][n]*Rlp[nyu][n] )*( Rlp[myu][n]*Rlp[nyu][n] )*((double)(n-1)*dr)*dr; }

					if ( modep[myu] == Ptotal || modep[nyu] == Ptotal ) { hmn = 0.0; }	// �ō������[�h�͖���
					else {
						hmn = deps2*((omega*omega*(PI*sqrt(PI))*(Dc*Dc*Dc)) / 8.0)*exp(-(dbmn*dbmn)*(Dc*Dc) / 4.0)*E_over; }

					if ( myu == nyu ) { H[myu][nyu] = 0.0; } else { H[myu][nyu] = hmn*dz; } // ������

					if ( y == Nl/2 || fout == 1 ) {
						if ( myu != nyu && myu > nyu ) {
//						if ( myu != nyu ) {
							fprintf ( fq, "%d,%d,%d,%d,%f,%d,%d,%d,%d,%f,%d,%f,%e,%e\n",	myu, modem [myu], model [myu], modep [myu], Mbeta [myu],
								nyu, modem [nyu], model [nyu], modep [nyu], Mbeta [nyu], abs(modep[myu]-modep[nyu]), fabs(dbmn), E_over, hmn ); }}
					E_over = 0.0; }}		
			}
			}
			// �}�C�N���x���f�B���O
			if ( wo == 2) {	dbmn = 0.0; hmn =0.0;
			for ( myu = 0; myu < NLP; myu++ ) {
				for ( nyu = 0; nyu < NLP; nyu++ ) {
					dbmn = Mbeta [myu] - Mbeta [nyu];

					if ( modep[myu] == Ptotal || modep[nyu] == Ptotal ) { hmn = 0.0; }	// �ō������[�h�͖���
					else if ( abs (modep [myu] - modep [nyu]) == 1 ) {
						if ( modep[myu] > modep [nyu] ) { mm = modep[myu] ; } else { mm = modep[nyu]; }
						hmn = (1.0/8.0)*(n1*k*A)*(n1*k*A)*pow( ((double)mm/(double)Ptotal), 4.0/(g+2.0) )
							*sigma2*sqrt(PI)*Db*exp(-(dbmn*dbmn)*(Db*Db) / 4.0); }
					else { hmn = 0.0; }

					if ( myu == nyu ) { H[myu][nyu] = 0.0; } else { H[myu][nyu] = hmn*dz; } // ������
					if ( y == Nl/2 || fout == 1 ) {
						if ( myu != nyu && myu > nyu ) {
//						if ( myu != nyu ) {
							fprintf ( fq, "%d,%d,%d,%d,%f,%d,%d, %d,%d,%f,%d,%f,%e\n", myu, modem [myu], model [myu], modep [myu], Mbeta [myu],
								nyu, modem [nyu], model [nyu], modep [nyu], Mbeta [nyu], abs(modep[myu]-modep[nyu]), fabs(dbmn), hmn ); }}
					 hmn = 0.0; }}
			}	//! H�s��̎Z�o�I��
			
			// ������		
			hmn = 0.0; myu = 0; nyu = 0;
			// �Ίi�v�f
			for ( myu = 0; myu < NLP; myu++ ) { Hrowsum = 0.0;
			for ( nyu = 0; nyu < NLP; nyu++ ) {
				if ( modem[nyu] == 0 ) { Hrowsum = Hrowsum + H[myu][nyu]; } else { Hrowsum = Hrowsum + 2.0*H[myu][nyu]; }}
			H[myu][myu] = 1.0 - ( 2.0*alpha[myu]*dz + Hrowsum ); }
			// ��Ίi�v�f
			for ( myu = 0; myu < NLP; myu++ ) {
			for ( nyu = 0; nyu < NLP; nyu++ ) {
				if ( myu != nyu ) { if ( modem[myu] == 0 ) { H[myu][nyu] = 1.0*H[myu][nyu]; } else { H[myu][nyu] = 2.0*H[myu][nyu]; }}}}
			// ��������̊m�F
			Hmmmin = H[0][0];
			for ( m = 1; m < NLP; m++ ) { if ( Hmmmin > H[m][m] ) { Hmmmin = H[m][m]; } }
			if ( Hmmmin < 0) { printf ( "Too large ��z value! Change the value appropriately!\n" ); exit ( EXIT_FAILURE ); }}
			/* H�s��̏o��*/
			if ( y == Nl/2 && fout == 1 ) {			// ������ҏW����ׂ���������܂���
				for ( myu = 0; myu < NLP; myu++ ) {	
				for ( nyu = 0; nyu < NLP; nyu++ ) { fprintf ( fr, "%f,", H[myu][nyu] ); } fprintf ( fr, "\n" ); }}

			//! ��U�����ݒ�iA+�s��̎Z�o�j
			int minput=0;
			// OFL condition
			if ( launch == 0 ) {
				for ( m = 0; m < NLP; m++ ) {
					if ( matdis == 0 ) { if (modem[m] == 0) { Amplu[minput][m][0] = 100.0; } else { Amplu[minput][m][0] = 200.0; }} // �k�ސ��̍l��
					if ( matdis == 1 ) { 
						if (modem[m] == 0) { Amplu[minput][m][0] = 100.0*OSmat[ (int)((lamda -lpmin)/dlp) ][0]; }
						else { Amplu[minput][m][0] = 200.0*OSmat[ (int)((lamda -lpmin)/dlp) ][0]; }}}}
			
			// RML condition
			if ( launch == 1 ) {
#pragma omp parallel			
				{
#pragma omp for
				for ( m = 0; m < NLP; m++ ) { Amev = Amod = 0.0;
					for ( i = 0; i < Nxy; i++ ) { xx = ( r0 - ( (double)Nxy*dx / 2.0 ) ) + (double)i*dx;
						for ( j = 0; j < Nxy; j++ ) { yy = - ( (double)Nxy*dy / 2.0 ) + (double)j*dy;
						rr = sqrt ( xx*xx + yy*yy ); nr = (int) (rr / dr) 	; // �؂�̂āH
						if ( nr == 0 ) { Rxy =  Rlp [m][0] + ( Rlp [m][1] - Rlp [m][0]) * (rr /dr); ;}
//						else { Rxy = Rlp [m][nr] + ( Rlp [m][nr+1] - Rlp [m][nr]) * ((rr - (double)nr*dr) /dr); } // interpolation
						else if ( nr <= N ) { Rxy = Rlp [m][nr] + ( Rlp [m][nr+1] - Rlp [m][nr]) * ((rr - (double)nr*dr) /dr); } // core (interpolation)
						else if ( nr > N ) { Rxy = Rlp [m][N]*( bessk (m, w*rr) / bessk (m, w*A) ); } // cladding
						Emev = Rxy*cos( (double)(modem[m])*atan2(yy, xx) );
						Emod = Rxy*sin( (double)(modem[m])*atan2(yy, xx) );
						Ein = exp ( - ((xx - r0)*(xx - r0) + yy*yy) / (w0*w0) ); // input field
						Amev = Amev + Emev*Ein*dx*dy; Amod = Amod + Emod*Ein*dx*dy; }}
						Amev = Amev*Amev / (((2.0*omega*Mu0)/Mbeta[m])*(PI*w0*w0/2.0));
						Amod =  Amod*Amod / (((2.0*omega*Mu0)/Mbeta[m])*(PI*w0*w0/2.0));
					if ( matdis == 0 ) { Amplu[minput][m][0] = 100.0* (Amod + Amev); }
					if ( matdis == 1 ) { Amplu[minput][m][0] = 100.0*( Amod + Amev )*OSmat[ (int)((lamda -lpmin)/dlp) ][0]; }
					//printf ("Amplu [%d][0] =%f\n", modem[m],Amplu [m][0] );
					}}}

			if (launch == 2) {
				FILE   *fp5, *fp8;			
				char   trash[65536];
				int    min, lin, i, j, Nxy, Nvin, couple, nrr1, nrr2, OFFres, OFFrange;
				int    *modem, *modemin, *model, *modelin;
				double k, A, n0, dr, r0, dx, dy, aa, w, Rinxy, tauin, betain, Rinfin, eigin;
				double Em_evin, Em_odin, Em_ev, Em_od, Am_evev, Am_evod, Am_odev, Am_odod, Avin;
				double *Mbeta, *Mbetain,    **Rlp, **MPD2d, **Rinlp;
				double win, xx1, xx2, yy1, yy2, rr1, rr2;
				double betaoverk, cef_odod, cef_evod, cef_odev, cef_evev, Adash = Vradius;
				if ((fp5 = fopen("[VCSEL_intensity_profile].csv", "r")) != NULL) {		//		VCSEL ��LP���[�h��1�������x���z�t�@�C�����J��
					fgets(trash, 65536, fp5);
					char MPDfilename[256];
					char coupletype[128];
					//! for���[�v  [OFFSET]
					//for (int ii = 0; ii <= OFFrange / OFFres + 1; ii++) {
						fscanf(fp5, "%s", trash);		//!? fp5�Ɋւ���fscanf�͔g�����[�v�̊O�ɂ��邩���C���֐����ɓ���邩������
						double E2m_evin, E2m_ev;
						//! for���[�v  [MODE NUMBER]
						for (minput = 0; minput < Nvmode; minput++) {
							
							fscanf(fp5, "%d, %d, %lf, %lf, %lf, %lf, %lf,", &min, &lin, &tauin, &betain, &betaoverk, &Rinfin, &eigin);
							modemin[minput] = min;
							modelin[minput] = lin;
							Mbetain[minput] = betain;
							
							printf("%d\t%d\n", modemin[minput], modelin[minput]);
							if (couple == 0) { Adash = Avin; }		if (couple == 1) { Adash = A; }

							for (j = 0; j <= (int)(Adash / dr); j++) {
								fscanf(fp5, "%lf,", &Rinlp[minput][j]);				}
							fscanf(fp5, "%lf\n", &Rinlp[minput][j]); 
							for (m = 0; m < NLP; m++) {
								Am_evev = 0.0;		Am_evod = 0.0;
								Am_odev = 0.0;		Am_odod = 0.0;
								E2m_evin = 0.0;		E2m_ev = 0.0;
								Rinxy = 0.0;

								for (i = 0; i < Nxy; i++) {
									xx1 = (Vradius / Avin) * ((-(double)Nxy * dx / 2.0) + ((double)i * dx));
									if (launch == 2) {
										xx2 = (r0 - (double)Nxy * dx / 2.0) + ((double)i * dx);
									}
									for (j = 0; j < Nxy; j++) {
										yy1 = (Vradius / Avin) * (-(double)Nxy * dx / 2.0) + ((double)j* dy);
										yy2 = (-(double)Nxy * dx / 2.0) + ((double)j * dy);
										rr1 = sqrt(xx1 * xx1 + yy1 * yy1);
										rr2 = sqrt(xx2 * xx2 + yy2 * yy2);
										nrr1 = (int)(rr1 / dr);
										nrr2 = (int)(rr2 / dr);
										//   �����͓��˂��郌�[�U�i�t�@�C�o�j�́Cxy���ʏ��2�������x���z���v�Z    
										if (nrr1 == 0) {		// �ȉ��̔z��Rinxy�́C���������̔z��Rlp���g���̂œǂݍ��݂�HDD�����葬��
											Rinxy = Rinlp[minput][0] + (Rinlp[minput][1] - Rinlp[minput][0]) * (rr1 / dr);
										}
										else if (nrr1 <= Nvin) {
											Rinxy = Rinlp[minput][nrr1] + (Rinlp[minput][nrr1 + 1] - Rinlp[minput][nrr1]) * ((rr1 - (double)nrr1 * dr) / dr);
										}
										else if (nrr1 > Nvin) {
											//printf("%d, %lf, %d, %lf\n", minput, Mbetain[minput], Nvin, Avin);
											win = aa * sqrt(Mbetain[minput] * Mbetain[minput] - k * k * n0 * n0);
											Rinxy = Rinlp[minput][Nvin] * (bessk(modemin[minput], win * rr1) / bessk(modemin[minput], win * Avin));
										}
										//   �����͎󂯎葤�̃t�@�C�o�ɂ�����Cxy���ʏ��2�������x���z���v�Z				/// �ȉ��̔z��Rinxy�́C���������̔z��Rlp���g���̂œǂݍ��݂�HDD�����葬��
										if (nrr2 == 0) {
											Rxy = Rlp[m][0] + (Rlp[m][1] - Rlp[m][0]) * (rr2 / dr);
										}
										else if (nrr2 <= N) {
											Rxy = Rlp[m][nrr2] + (Rlp[m][nrr2 + 1] - Rlp[m][nrr2]) * ((rr2 - (double)nrr2 * dr) / dr);
										}
										else if (nrr2 > N) {
											w = aa * sqrt(Mbeta[m] * Mbeta[m] - k * k * n0 * n0);
											Rxy = Rlp[m][N] * (bessk(modem[m], w * rr2) / bessk(modem[m], w * A));
										}

										//   �d�Ȃ�ϕ��ɗp����d�ꕪ�z�̎Z�o�E�d�Ȃ�ϕ��̎��s   
										Em_evin = Rinxy * cos((double)(modemin[minput]) * atan2(yy1, xx1));
										Em_odin = Rinxy * sin((double)(modemin[minput]) * atan2(yy1, xx1));
										Em_ev = Rxy * cos((double)(modem[m]) * atan2(yy2, xx2));
										Em_od = Rxy * sin((double)(modem[m]) * atan2(yy2, xx2));
										if (i < 1 && j < 1) {
											if (modem[m] == 0) {
												//printf("%lf, %lf, %lf, %lf\t\t", Em_evin, Em_odin, Em_ev, Em_od); 
												printf("%d ", m);
											}
											else { printf("%d ", m); }
										}
										E2m_evin = E2m_evin + Em_evin * Em_evin * dx * dy;
										E2m_ev = E2m_ev + Em_ev * Em_ev * dx * dy;
										Am_evev = Am_evev + Em_evin * Em_ev * dx * dy;
										Am_evod = Am_evod + Em_evin * Em_od * dx * dy;
										Am_odev = Am_odev + Em_odin * Em_ev * dx * dy;
										Am_odod = Am_odod + Em_odin * Em_od * dx * dy;
									}
								}
								cef_evev = Am_evev * Am_evev / (E2m_evin * E2m_ev);
								cef_evod = Am_evod * Am_evod / (E2m_evin * E2m_ev);
								cef_odev = Am_odev * Am_odev / (E2m_evin * E2m_ev);
								cef_odod = Am_odod * Am_odod / (E2m_evin * E2m_ev);

								if (modem[m] == 0 && modemin[minput] == 0) {
									MPD2d[m][minput] = 100.0 * (cef_evev + cef_evod + cef_odev + cef_odod);
								}
								if (modem[m] != 0 && modemin[minput] == 0) {
									MPD2d[m][minput] = 100.0 * (cef_evev + cef_evod + cef_odev + cef_odod) / 2;
								}
								if (modem[m] == 0 && modemin[minput] != 0) {
									MPD2d[m][minput] = 100.0 * (cef_evev + cef_evod + cef_odev + cef_odod) / 2;
								}
								if (modem[m] != 0 && modemin[minput] != 0) {
									MPD2d[m][minput] = 100.0 * (cef_evev + cef_evod + cef_odev + cef_odod) / 2;
								}
							}
							printf("\n");
						}
						//! Almup�Ɍv�Z���ʂ���
						for (minput = 0; minput < Nvmode; minput++) {
							for (m = 0; m < NLP; m++) {
								Amplu[minput][m][0] = MPD2d[m][minput];	}
						}
						printf("\n");
						for (minput = 0; minput < Nvmode; minput++) { printf("%d\t%lf\n", modemin[minput], Mbetain[minput]); }
						printf("\n");
						for (m = 0; m < NLP; m++) { printf("%d\t%lf\n", modem[m], Mbeta[m]); }
						if (launch == 2) { sprintf(MPDfilename, "MPD_%s\\MPD2d_%s_%dum_Offset_single(not_dependence_on_x).csv", coupletype, coupletype, r0); }
						if ((fp8 = fopen(MPDfilename, "w")) != NULL) {
							fprintf(fp8, "Mode number for receiver-side fiber,");
							if (couple == 0) { fprintf(fp8, "m, l, "); }
							for (j = 0; j < Nvmode; j++) { fprintf(fp8, "%d,", j); }
							fprintf(fp8, "\n");
							for (i = 0; i < NLP; i++) {
								fprintf(fp8, "%d,", i);
								if (couple == 0) { fprintf(fp8, "%d, %d, ", modem[i], model[i]); }
								for (j = 0; j < Nvmode; j++) {
									// printf("%e, ", MPD2d[i][j]);
									fprintf(fp8, "%e,", MPD2d[i][j]);
								}
								fprintf(fp8, "\n");
							}
							fclose(fp8);
						}
						if (launch == 2) break;
					//}
				}

				/*for (n = 0;n < Nom; n++) {	
				{
					for ( m = 0; m < NLP; m++ ) { Amev = Amod = 0.0;
						for ( i = 0; i < Nxy; i++ ) { xx = ( r0 - ( (double)Nxy*dx / 2.0 ) ) + (double)i*dx;
							for ( j = 0; j < Nxy; j++ ) { yy = - ( (double)Nxy*dy / 2.0 ) + (double)j*dy;
							rr = sqrt ( xx*xx + yy*yy ); nr = (int) (rr / dr) 	; // �؂�̂āH
							if ( nr == 0 ) { Rxy =  Rlp [m][0] + ( Rlp [m][1] - Rlp [m][0]) * (rr /dr); ;}
	//						else { Rxy = Rlp [m][nr] + ( Rlp [m][nr+1] - Rlp [m][nr]) * ((rr - (double)nr*dr) /dr); } // interpolation
							else if ( nr <= N ) { Rxy = Rlp [m][nr] + ( Rlp [m][nr+1] - Rlp [m][nr]) * ((rr - (double)nr*dr) /dr); } // core (interpolation)
							else if ( nr > N ) { Rxy = Rlp [m][N]*( bessk (m, w*rr) / bessk (m, w*A) ); } // cladding
							Emev = Rxy*cos( (double)(modem[m])*atan2(yy, xx) );
							Emod = Rxy*sin( (double)(modem[m])*atan2(yy, xx) );
							Ein = exp ( - ((xx - r0)*(xx - r0) + yy*yy) / (w0*w0) ); // input field
							Amev = Amev + Emev*Ein*dx*dy; Amod = Amod + Emod*Ein*dx*dy; }}
							Amev = Amev*Amev / (((2.0*omega*Mu0)/Mbeta[m])*(PI*w0*w0/2.0));
							Amod =  Amod*Amod / (((2.0*omega*Mu0)/Mbeta[m])*(PI*w0*w0/2.0));
						if ( matdis == 0 ) { Amplu[minput][m][0] = 100.0* (Amod + Amev); }
						if ( matdis == 1 ) { Amplu[minput][m][0] = 100.0*( Amod + Amev )*OSin[ (int)((lamda -lpmin)/dlp) ]; }
						//printf ("Amplu [%d][0] =%f\n", modem[m],Amplu [m][0] );
						}}}*/
				
			}
			//TODO End of lanching condition setting for measuring minEMBc   //////////////////////////////////////////////
			/*if ( launch == 2 ) { E2total = 0.0;
				for ( m = 0; m < NLP; m++ ) { Amev = Amod = 0.0;
					for ( i = 0; i <= Nxy; i++ ) { xx = - ( (double)Nxy*dx / 2.0 ) + (double)i*dx; // ���S��U
						for ( j = 0; j <= Nxy; j++ ) { yy = - ( (double)Nxy*dy / 2.0 ) + (double)j*dy;
						rr = sqrt ( xx*xx + yy*yy ); nr = (int) (rr / dr) 	; // �؂�̂āH
						if ( nr == 0 ) { Rxy =  Rlp [m][0] + ( Rlp [m][1] - Rlp [m][0]) * (rr / dr); ;}
						else if ( nr <= N ) { Rxy = Rlp [m][nr] + ( Rlp [m][nr+1] - Rlp [m][nr]) * ((rr - (double)nr*dr) /dr); } // core (interpolation)
						else if ( nr > N ) { w = aa*sqrt( Mbeta[m]*Mbeta[m] - k*k*n0*n0 ); Rxy = Rlp [m][N]*( bessk (m, w*rr) / bessk (m, w*A) ); } // cladding
						Emev = Rxy*cos( (double)(modem[m])*atan2(yy, xx) );
						Emod = Rxy*sin( (double)(modem[m])*atan2(yy, xx) );
						Ein = sqrt (Pxyin[i][j]); // input field
						if (m == 0) { E2total = E2total + Ein*Ein*dx*dy; }
						Amev = Amev + Emev*Ein*dx*dy; Amod = Amod + Emod*Ein*dx*dy; }}
						Amev = Amev*Amev / ((2.0*omega*Mu0)/(PI*Mbeta[m]));
						Amod =  Amod*Amod / ((2.0*omega*Mu0)/(PI*Mbeta[m]));
					    Aplu[m][0] = 100.0* (Amod + Amev) / E2total;
					//printf ("Aplu[%d][0]=%f Amev=%f Amod=%f E2total=%f\n", modem[m],Aplu [m][0],Amev,Amod,E2total );
					} free_dmatrix ( Pxyin, 0, Nxy, 0, Nxy ); }*/

			// This may be not adequate. for ( m = 0; m < NLP; m++ ) { A00 = A00 + Aplu[m][0]; } Aw00 = Aw00 + A00;
			// Correction for the highest mode group ( elimination of power )
			for (minput = 0;minput < Nvmode;minput++) {
				for ( m = 0; m < NLP; m++ ) { if ( modep[m] == Ptotal ) { Amplu[minput][m][0] = 0.0; }}}
			// Total power of input impulse
			for(minput=0;minput<Nvmode;minput++){
				for ( m = 0; m < NLP; m++ ) { A00 = A00 + Amplu[minput][m][0]; } Aw00 = Aw00 + A00;}

			fprintf ( fp2, "Amplitude of input impulse,A00,%e,\n", A00 );
			fprintf ( fp2, "Cumulative amplitude of input impulse,Aw00,%e,\n\n", Aw00 );
			free_dmatrix ( Rlp, 0, 2*Nwkb, 0, N-1 );
			
			/* �o�̓t�@�C�� */
			if ( y == Nl/2 || fout == 1 ) {
				// ���Ԕg�` P
				fprintf ( fr2, "nstd=%d,", nstd );
				for ( n = 0; n <= Lmax; n++ ) { fprintf ( fr2, "%f,", (double)n*Tv );} fprintf ( fr2, "pct.\n" );
				fprintf ( fr2, "%f,%f\n", 0.0, A00 );
				// ���[�h�p���[���z Pm
				fprintf ( fr3, "myu," ); for ( m = 0; m < NLP; m++ ) { fprintf ( fr3, "%d,", m ); } fprintf ( fr3, "\n" );
				fprintf ( fr3, "LP," ); for ( m = 0; m < NLP; m++ ) { fprintf ( fr3, "LP%d_%d,", modem[m], model[m] ); } fprintf ( fr3, "\n" );
				fprintf ( fr3, "m," ); for ( m = 0; m < NLP; m++ ) { fprintf ( fr3, "%d,", modem[m] ); } fprintf ( fr3, "\n" );
				fprintf ( fr3, "l," ); for ( m = 0; m < NLP; m++ ) { fprintf ( fr3, "%d,", model[m] ); } fprintf ( fr3, "\n" );
				fprintf ( fr3, "p," ); for ( m = 0; m < NLP; m++ ) { fprintf ( fr3, "%d,", modep[m] ); } fprintf ( fr3, "\n" );
				fprintf ( fr3, "%f,", 0.0 );
				for ( m = 0; m < NLP; m++ ) { for ( n = 0; n <= Lmax; n++ ) { Pm[m] = Pm[m] + Amplu[minput][m][n]; }
				fprintf ( fr3, "%f,", Pm[m] ); } fprintf ( fr3,"\n"); 
				// ���[�h�Q�p���[���z Pg
				fprintf ( fr9, "," ); for ( j = 1; j <= Ptotal; j++ ) { fprintf ( fr9, "%d,", j ); } fprintf ( fr9, "\n" ); // ���[�h�Q�ԍ�					
				fprintf ( fr9, "," ); for ( i = 1; i <= Ptotal; i++ ) { count =0; 
				for ( myu = 0; myu < NLP; myu++ ) { if ( modep [myu] == i ) { count = count +1; }} // �k�ސ��i���[�h���j 
				fprintf ( fr9, "%d,", count ); pdeg[i-1] =count; } fprintf ( fr9,"\n");
				fprintf ( fr9, "%f,", 0.0 ); 
				for ( m = 0; m < NLP; m++ ) { Pg[modep[m]-1] = Pg[modep[m]-1] + Pm[m]; } // ���[�h�Q�p���[
				for ( j = 0; j < Ptotal; j++ ) { fprintf ( fr9, "%f,", Pg[j] / (double)pdeg[j] ); Pg[j] = 0.0; }  fprintf ( fr9,"\n");			
			}

			////! /////        CONTINUE��        //////
			////printf("%lf\t", OSin[(int)((lamda - lpmin) / dlp)]);
			////	if (OSin[(int)((lamda - lpmin) / dlp)] == 0.0) {
			////	//contct += 1;		printf("continue count is %d\n", contct);
			////	continue;	}

			//! ////////////////////////////////////////////////////////
			//! //////////           5. ���g�`�����             //////////
			//! ////////////////////////////////////////////////////////
			/****************************************************************************************/
			printf("start!\n");
			if(launch ==0 || launch==1){ Nvmode=1;}
			minput = 0;
			for ( i = 1; i <= Nz; i++ ) // z(i-1) ~ zi
			{
				/* ���Βx���X�e�b�v���̍ő�l�Z�o */
				Li =0;
				Li = (int) ( (( taumax - taumin )*( (double)i*dz ) / Tv) + 0.5 );				

				/* ���Βx���X�e�b�v���̃��[�h���z�Z�o */
				fprintf(fdebug, "kim[m], ");
				for ( m = 0; m < NLP; m++ ) {
					kim[m] = (int) ( (( Mtau[m] - taumin )*( (double)i*dz ) / Tv ) + 0.5 )
					              - (int) ( (( Mtau[m] - taumin )*( (double)(i-1)*dz ) / Tv ) + 0.5 );
					fprintf(fdebug, "%d, ", kim[m]);}	
				fprintf(fdebug, "\n");

#pragma omp parallel
{
#pragma omp for
				/* �^�C���V�t�g���Z */
				for (minput=0; minput<Nvmode; minput++){
					for ( m = 0; m < NLP; m++ ) {
						for ( n = 0; n < kim[m]; n++ ) { Amin[minput][m][n] = 0.0; }
						for ( n = kim[m]; n <= Li; n++ ) { Amin[minput][m][n] = Amplu[minput][m][n - kim[m]]; }}}
#pragma omp for
				/* �J�b�v�����O���Z */
				for (minput = 0; minput < Nvmode; minput++) {					
					for ( m = 0; m < NLP; m++ ) {
						for ( n = 0; n <= Li; n++ ) { me=0.0;
						for ( l = 0; l < NLP; l++ ) { me = me + H[m][l]*Amin[minput][l][n]; }
						Amplu[minput][m][n] = me; }}}
}

				// �@ �g�����U���l������ꍇ
				/* ���Ԕg�`�̏d�ˍ��킹 */
				minput = 0;
				if ( matdis == 1 && i == Nz ) {
					// �e�g�������̃C���p���X�����i����Ԃ͔�l���j
					if ( y == Nl/2 || fout == 1) { fprintf ( fr2, "%f,", (double)i*dz ); }
					for ( n = 0; n <= Li; n++ ) { ap = 0.0;
						for ( m = 0; m < NLP; m++ ) { ap = ap + Amplu[minput][m][n]; }
						if ( y == Nl/2 || fout == 1) { 
							fprintf ( fr2, "%f,", ap ); 
							if(n==0){ fprintf(fdebug, "ap,"); }
							fprintf ( fdebug, "%f,", ap );}}
					fprintf(fdebug, "\n");
					if ( y == Nl/2 || fout == 1) { fprintf ( fr2,"\n" ); }
					//! �ŏ��Q�x���̔g���ˑ������l�������������킹�i����Ԃ��l���j
					if ( y == 0 || nstd == 0 ) { 
						for ( n = 0; n <= Li; n++ ) { 
							for ( m = 0; m < NLP; m++ ) { 
								P[n] = P[n] + Amplu[minput][m][n]; }}; }
					if ( y != 0 && nstd > 0 ) { 
						for ( n = 0; n <= Li; n++ ) { 
							for ( m = 0; m < NLP; m++ ) { 
								P[n+nstd] = P[n+nstd] + Amplu[minput][m][n]; }}; }
					if ( y != 0 && nstd < 0 ) {
						if ( nstd < nstdmin ) {
							for ( n = 0; n <= Pnum[y-1] ; n++ ) { 
								P[(Pnum[y-1]-n)+(nstdmin-nstd)] = P[(Pnum[y-1]-n)]; }
							for ( n = 0; n < nstdmin-nstd; n++ ) { P[n] = 0.0; }
							for ( n = 0; n <= Li; n++ ) { for ( m = 0; m < NLP; m++ ) { 
								P[n] = P[n] + Amplu[minput][m][n]; }}}
						else {
							for ( n = 0; n <= Li; n++ ) { for ( m = 0; m < NLP; m++ ) {
								P[n+(nstd-nstdmin)] = P[n+(nstd-nstdmin)] + Amplu[minput][m][n]; }}; }};
					// for ( n = 0; n <= Pnum[y]; n++ ) { fprintf ( fs3,"%f,", P[n] ); } fprintf ( fs3,"\n" );
				}

#if 0
				// �A �g�����U���l�����Ȃ��ꍇ	
				/* �v�Z���ʂ̏o�� */	
				if ( matdis == 0 && (i%Nzout == 0 || i == Nz) ) {
					/* ���[�h�p���[���z Pm */
					if ( fout == 1 || y == Nl/2 ) { fprintf ( fr3, "%f,", (double)i*dz ); }
					for ( m = 0; m < NLP; m++ ) {
						Pm[m] = 0.0;
						for ( n = 0; n <= Li; n++ ) { Pm[m] = Pm[m] + Amplu[minput][m][n]; }
						if ( fout == 1 || y == Nl/2 ) { fprintf ( fr3, "%f,", Pm[m] ); }					
						Pg[modep[m]-1] = Pg[modep[m]-1] + Pm[m]; }
					if ( fout == 1 || y == Nl/2 ) { fprintf ( fr3,"\n"); }
					/* ���[�h�Q�p���[���z Pg */
					if ( fout == 1 || y == Nl/2 ) { fprintf ( fr9, "%f,", (double)i*dz );						
					for ( j = 0; j < Ptotal; j++ ) { fprintf ( fr9, "%f,", Pg[j] / (double)pdeg[j] ); Pg[j] =0.0;}  fprintf ( fr9,"\n"); }
					/* �C���p���X���� P */
					if ( fout == 1 || y == Nl/2) { fprintf ( fr2, "%f,", (double)i*dz ); }
					for ( n = 0; n <= Li; n++ ) {
						P[n] = 0.0;
						for ( m = 0; m < NLP; m++ ) { P[n] = P[n] + Amplu[minput][m][n]; }														
						if ( fout == 1 || y == Nl/2 ) { fprintf ( fr2, "%f,", P[n] ); } }
					if ( fout == 1 || y == Nl/2 ) { fprintf ( fr2,"\n" ); }										
					/* �o�͔g�` */						 
					Pout = drealvector ( 0, Ti+Li ); init_realvector ( Pout, 0, Ti+Li );
					Pout = convolution ( Ti, Pin, Li, P );
					if ( fout == 1 || y == Nl/2) {
						for ( n = 0; n < Ti+Li; n++ ) { fprintf ( fr6, "%f,", Pout[n] );} fprintf ( fr6,"\n" );}
					/* ���g������ M */
					if ( fout == 1 || y == Nl/2 ) { fprintf ( fr4, "%f,", (double)i*dz ); }
					for ( j = 0; j <= Nf; j++ ) {
						Hw = ReHw = ImHw = 0.0;
						for ( n = 0; n <= Li; n++ ) {
							ReHw = ReHw + P[n]*cos ( (2.0*PI*(fmin+(double)j*df))*(double)n*Tv );
							ImHw = ImHw - P[n]*sin ( (2.0*PI*(fmin+(double)j*df))*(double)n*Tv ); }
						M[j] = sqrt ( ReHw*ReHw + ImHw*ImHw ) / A00; if ( y == Nl/2 || fout == 1 ) { fprintf ( fr4, "%f,", -10.0*log10(1.0 / M[j]) ); }}
					/* -3dB�ш敝 bw */
					bw = tav = Ptot = rms = 0.0;
					for ( j = 0; j <= Nf-1; j++ )  { if ( M[j] >0.5 && M[j+1] < 0.5 ) { bw = (fmin*1.0e3) + (df *1.0e3)*(j + (M[j] - 0.5) / (M[j] - M[j+1])); break; } }
					if ( fout == 1 || y == Nl/2 ) { fprintf ( fr4, "%f,", bw );}		
					/* �C���p���X����RMS�� rms */
					for ( n = 0; n <= Li; n++ ) { tav = tav + ( taumin*(double)(i-1)*dz + (double)n*Tv)*P[n]*Tv; Ptot = Ptot + P[n]*Tv; } tav = tav / Ptot;
					for ( n = 0; n <= Li; n++ ) { rms = rms + ( taumin*(double)(i-1)*dz + (double)n*Tv)*( taumin*(double)(i-1)*dz + (double)n*Tv)*( P[n] / Ptot )*Tv; }
					rms = sqrt ( rms - (tav*tav) );	

					if ( fout == 1 || y == Nl/2 ) { fprintf ( fr4, "%f,%f", rms, tav - ( taumin*(double)(i-1)*dz) ); fprintf ( fr4, "\n" ); }

					fprintf ( fs, "%f,%d,%d,%f,%f,%f,%f,%f,%f,%f\n", g, NLP, NLP0, taumin*(double)(i-1)*dz*1.0e-3, (double)i*dz, rms, bw, tav, Ptot, Aw00 );
					printf ( "length:%f m ", (double)i*dz ); printf ( "-3db bandwidth:%f GHz ", bw ); printf ( "pulse broadening:%f ps\n", rms );
				}	
#endif
			}
			/****************************************************************************************/

			 if ( nstd < nstdmin ) { nstdmin = nstd; }
			 if ( nstd > nstdmax ) { nstdmax = nstd; }

			 free_dmatrix ( H, 0, NLP-1, 0, NLP-1 );
			 free_d3tensor ( Amin, 0, Nvmode, 0, NLP-1, 0, Lmax );
			 free_d3tensor ( Amplu, 0, Nvmode, 0, NLP-1, 0, Lmax );
			 free_drealvector ( alpha, 0 );
			 free_dintvector ( kim, 0 );
			 free_drealvector ( Pm, 0 );
			 free_drealvector ( Pg, 0 );
			 free_dintvector ( modem, 0 );
			 free_dintvector ( model, 0 );
			 free_dintvector ( modep, 0 );
			 free_dintvector ( pdeg, 0 );
			 free_drealvector ( Mbeta, 0 );
			 free_drealvector ( Mtau, 0 );
			 printf ( "A00=%f\n", A00 );
			 printf ( "End\n\n");			 
} //! �g�����[�v�I��
printf ( "Aw00=%f\n", Aw00 );

//! ///////////////////////////////////////////////////////
//! //////////           6. ���ʂ̏o��             //////////
//! //////////////////////////////////////////////////////

#if 1
/* �v�Z���ʂ̏o�� */	
if ( matdis == 1 ) {
	/* �C���p���X�����g�`�i�����X�y�N�g���l���j Pw */
	for ( n = 0; n <= Pnum[Nl]; n++ ) { fprintf ( fs3,"%f,", P[n] ); } fprintf ( fs3,"\n" );
	/*�@�o�͔g�` */
	Pout = drealvector ( 0, Ti+Pnum[Nl] ); init_realvector ( Pout, 0, Ti+Pnum[Nl] );	
	Pout = convolution ( Ti, Pin, Pnum[Nl], P );
	if ( fout == 1 || y == Nl/2) {
		for ( n = 0; n < Ti+Pnum[Nl]; n++ ) { fprintf ( fs4, "%f,", Pout[n] );} fprintf ( fs4,"\n" );}	
	/* ���g������ M */	
	fprintf ( fr4, "%f,", (double)Nz*dz );
	for ( j = 0; j <= Nf; j++ ) {
		ReHw = ImHw = 0.0;
		for ( n = 0; n <= Pnum[Nl]; n++ ) {
			ReHw = ReHw + P[n]*cos ( (2.0*PI*(fmin+(double)j*df))*(double)n*Tv );
			ImHw = ImHw - P[n]*sin ( (2.0*PI*(fmin+(double)j*df))*(double)n*Tv ); }
		M[j] = sqrt ( ReHw*ReHw + ImHw*ImHw ) / Aw00; fprintf ( fr4, "%f,", -10.0*log10(1.0 / M[j]) ); }

	/* -3dB�ш敝 bw */
	bw = tav = Ptot = rms = 0.0;	
	for ( j = 0; j <= Nf-1; j++ )  { if ( M[j] >0.5 && M[j+1] < 0.5 ) { bw = (fmin*1.0e3) + (df *1.0e3)*(j + (M[j] - 0.5) / (M[j] - M[j+1])); break; }}
	fprintf ( fr4, "%f,", bw );

	/* �C���p���X����RMS�� rms */
	for ( n = 0; n <= Pnum[Nl]; n++ ) {	tav = tav + ( taumin*(double)(i-1)*dz + (double)n*Tv)*P[n]*Tv; Ptot = Ptot + P[n]*Tv; } tav = tav / Ptot;
	for ( n = 0; n <= Pnum[Nl]; n++ ) {	rms = rms + ( taumin*(double)(i-1)*dz + (double)n*Tv)*( taumin*(double)(i-1)*dz + (double)n*Tv )*( P[n] / Ptot )*Tv; }						
	rms = sqrt ( rms - (tav*tav) );
	fprintf ( fr4, "%f,", rms); fprintf ( fr4, "\n" );
					
	fprintf ( fs, "%f,%f,%f,%f,%f,%f,%f,", g, (double)Nz*dz, rms, bw, tav, Ptot, Aw00 );
	printf ( "length:%f m\n", (double)Nz*dz ); printf ( "-3db bandwidth:%f GHz\n", bw ); printf ( "rms width:%f ps\n", rms ); }

#endif


/* �X�y�b�N���R���g���X�g�̌v�Z */	
#if 0
if ( scc == 1 ) {
	/* �����X�y�N�g�����ȑ��֊֐� Cp */
	Cp = drealvector ( 0, Nf ); init_realvector ( Cp, 0, Nf );
	for ( j = 0; j <= Nf; j++ ) {
		for ( i = 0; i <= Nfp; i++ ) {
			int jj = (int) ( (double)j*df / dfp );
			Cp[j] = Cp[j] + OSmat[Nfp-i][0]*OSmat[ Nfp-i+jj ][0]; }
		Cpsum = Cpsum + df*Cp[j]; }


// �K�E�X�^�X�y�N�g���`��֐�
//		Cp[j] = Cp[j] + exp( -pow( (fpmin + dfp*double(i) - fp0) / fpgw, 2.0 ))
//								*exp( -pow( (fpmin + dfp*double(i) - df*double(j) - fp0 ) / fpgw, 2.0 ))*dfp; }

	/* �X�y�b�N���R���g���X�g spct */
	fprintf ( fr4, ","  );
	for ( j = 0; j <= Nf; j++ ) {	
		Cp[j] = Cp[j] / (2.0*Cpsum); fprintf ( fr4, "%f,", Cp[j] ); // Cp����֐��ł��邱�Ƃ��l��
		spct = spct + Cp[j]*M[j]*M[j]*df; } spct = sqrt ( 2.0*spct );
		fprintf ( fs, "%f\n", spct );	printf ( "speckle contrast:%f\n", spct );
		free_drealvector (Cp, 0); }
#endif

#if 1
////printf("ifdis=1\n");
free_dintvector ( Pnum, 0 );
free_drealvector ( P, 0 );	
free_drealvector ( Pin, 0 );
free_drealvector ( Pout, 0 );
free_drealvector ( M, 0 );
#endif

fclose ( fp2 );	// BW_setting.csv
fclose ( fp3 );	// BW_profile.csv
fclose ( fp4 );	// FEM_result.csv
fclose (fq);	// Hmn_result.csv
fclose ( fr );	// CPE_Hmatrix.csv
fclose ( fr2 );	// CPE_Impulse-responce.csv			
fclose ( fr3 );	// CPE_Mode-power-distribution.csv
fclose ( fr4 );	// CPE_Frequency-response.csv
fclose ( fr6 );	// CPE_Output-pulse-waveform.csv
fclose (fs);	// BW_result.csv
fclose (fs2);	// BW_Source-spectrum.csv
fclose (fs3);	// BW_Impulse-responce.csv
fclose (fs4);	// BW_Output-pulse-waveform.csv
fclose (fdebug);
system("pause");
return 0;
}

void selectOsciMode(int Nvmode, double Vradius) {	//!? m,l,tau,beta,betaoverk,Rinf,eig,InfProf�������ɂ��Ă��悢���s�K�v�Ȃ̂ŗv����
	
	//! �錾�Ə�����
	FILE* fptr, *fptrcol, *fw;
	char fppath[128], st[65536];
	char s1[128];
	int    rcnt, colcnt;
	double dr;
	double lf1 = 0.0;
	int    *m_guided, *l_guided, *modenum;
	double *tau_guided, *beta_guided, * betaoverk_guided, * Rinf_guided, * eig_guided;
	double **IntProf_guided;
	int    *m, *l;
	double *tau, *beta, *betaoverk, *Rinf, *eig;
	double **IntProf;	
	
	sprintf(fppath, "FILE_for_IO\\FEM_result_vcsel.csv");
	
	//! FEM_result_vcsel.csv�́i�񐔁E�j�s���̊m�F
	fptrcol = fopen(fppath, "r");
	printf("%s\n", fppath);
	fscanf(fptrcol, "%[^,], %[^,], %[^,], %[^,], %[^,], %[^,], %[^,],", &s1, &s1, &s1, &s1, &s1, &s1, &s1);
	// �񐔂��ȉ��ŃJ�E���g
	rcnt = 0;
	for (int i = 0 ; ; i++) {
		if (i == 0) { 
			fscanf(fptrcol, "%[^,],", &s1); 
			printf("%s, ", s1);		}
		else if (i == 1){
			fscanf(fptrcol, "%lf,", &dr);	//!? dr[um]
			printf("%lf, ", dr);	}
		else{ 
			fscanf(fptrcol, "%lf,", &lf1); 
			printf("%lf, ", lf1);	}		
		if (lf1 != Vradius) { rcnt += 1; }	//!? 3.5 um (VCSEL�R�A�N���b�h�E��)��ǂݍ��ނ܂Ōp��
		else { break; }						//!? �ǂݍ��݂������ubreak�v
		if (rcnt >= 10000) { break; system("pause"); }
	}
	printf("\n");

	// �s�����ȉ��ŃJ�E���g
	colcnt = -1;		//!? ���𒲐�
	while (fgets(st, 65536, fptrcol) != NULL) { colcnt++; }

	//! �s���E�񐔂̏o��
	printf("FEM_result_vcsel.csv\n");
	printf("Number of spatial step size for intensity profile: %d\n", rcnt);
	printf("Number of guided-mode in VCSEL cavity (LP): %d\n", colcnt);
	fclose(fptrcol);


	// �x�N�g���E�s�� �� �錾&������
	m_guided = dintvector(0, colcnt);				init_intvector(m_guided, 0, colcnt);
	l_guided = dintvector(0, colcnt);				init_intvector(l_guided, 0, colcnt);
	tau_guided = drealvector(0, colcnt);			init_realvector(tau_guided, 0, colcnt);
	beta_guided = drealvector(0, colcnt);			init_realvector(beta_guided, 0, colcnt);
	betaoverk_guided = drealvector(0, colcnt);		init_realvector(betaoverk_guided, 0, colcnt);
	Rinf_guided = drealvector(0, colcnt);			init_realvector(Rinf_guided, 0, colcnt);
	eig_guided = drealvector(0, colcnt);			init_realvector(eig_guided, 0, colcnt);
	IntProf_guided = dmatrix(0, colcnt, 0, rcnt);	init_dmatrix(IntProf_guided, 0, colcnt, 0, rcnt);
	// ���U���[�h�I���̍ۂ̉��Z�Ɏg�p
	modenum = dintvector(0, colcnt);		init_intvector(modenum, 0, colcnt);

	//!? //////    �ȉ����C���֐��Ɉړ����Ă��悢���s�K�v�Ȃ̂ŗv����    //////
	//!? m, l, tau, beta, betaoverk, Rinf, eig, InfProf �������ɂ���
	m = dintvector(0, Nvmode);				init_intvector(m, 0, Nvmode);
	l = dintvector(0, Nvmode);				init_intvector(l, 0, Nvmode);
	tau = drealvector(0, Nvmode);			init_realvector(tau, 0, Nvmode);
	beta = drealvector(0, Nvmode);			init_realvector(beta, 0, Nvmode);
	betaoverk = drealvector(0, Nvmode);		init_realvector(betaoverk, 0, Nvmode);
	Rinf = drealvector(0, Nvmode);			init_realvector(Rinf, 0, Nvmode);
	eig = drealvector(0, Nvmode);			init_realvector(eig, 0, Nvmode);
	IntProf = dmatrix(0, Nvmode, 0, rcnt);	init_dmatrix(IntProf, 0, Nvmode, 0, rcnt);
	
	if ((fptr = fopen(fppath, "r")) != NULL) {	
		fgets(st, 65536, fptr);
		for (int i = 0; i < colcnt; i++) {			
			fscanf(fptr, "%d, %d, %lf, %lf, %lf, %lf, %lf,", &m_guided[i], &l_guided[i], &tau_guided[i], 
				&beta_guided[i], &betaoverk_guided[i], &Rinf_guided[i], &eig_guided[i]);
			printf("%d, %d, %lf, %lf, %lf, %lf, %lf\n", m_guided[i], l_guided[i], tau_guided[i],
				beta_guided[i], betaoverk_guided[i], Rinf_guided[i], eig_guided[i]);
			for (int j = 0; j <= rcnt; j++) {
				fscanf(fptr, "%lf,", &IntProf_guided[i][j]);	}}
	}
	else { printf(" U cannot open the file !\n"); exit(EXIT_FAILURE); }
	
	fw = fopen("[VCSEL_intensity_profile].csv", "w");
	printf("m, l, tau, beta, ne, eig, R[N]\n");
	fprintf(fw, "m,l,tau,beta,ne,R_infinite,eig,r=0_um,");
	for (int i = 1; i <= rcnt; i++) {
		fprintf(fw, "%.3lf,", i*dr);	}//for (j = 1; j <= N; j++) { fprintf(fp4, "%lf,", (double)j * dr * 1.0e6); } fprintf(fp4, "\n");
	fprintf(fw, "\n");

	// ���U���[�h��I�肷�镔��
	for (int i = 0; i < colcnt;i++) {
		for (int j = 0;j < colcnt;j++) {
			if (beta_guided[i] < beta_guided[j]) {
				modenum[i] += 1;	}}}

	//! ���U���[�h��`���萔�̑傫�����ɕ��בւ��A�t�@�C��fw�ɏo��
	for (int i = 0; i < Nvmode;i++) {
		for (int j = 0;j < colcnt;j++) {
			if (modenum[j] == i) {
				m[i] = m_guided[j];
				l[i] = l_guided[j];
				tau[i] = tau_guided[j];
				beta[i] = beta_guided[j];
				betaoverk[i] = betaoverk_guided[j];
				Rinf[i] = Rinf_guided[j];
				eig[i] = eig_guided[j];
				fprintf(fw, "%d, %d, %lf, %lf, %lf, %lf, %lf, ", 
					m[i], l[i], tau[i], beta[i], betaoverk[i], Rinf[i], eig[i]);
				for (int k = 0; k <= rcnt; k++) {
					IntProf[i][k] = IntProf_guided[j][k];	
					fprintf(fw, "%lf, ", IntProf[i][k]);	}
				fprintf(fw, "\n");		}}
	}
	fclose(fw);
}

void inputFEM() {
	FILE   *fp, *fp2, *fp3, *fp4, *fp5, *fp6, *fr, *fpulse, *fopulse, *fmpd;
	char   fppath[128], fp2path[128], fp3path[128], fp4path[128], fp5path[128], fp6path[128], fp7path[128];
	char   fp8path[128], frpath[128], fpulsepath[128], fopulsepath[128], fmpdpath[128];
	          fppath[0]   =	 fp2path[0] = fp3path[0] = fp4path[0] = fp5path[0] = fp6path[0] = fp7path[0] =
	          fp8path[0]  =  frpath[0] = fpulsepath[0] = fopulsepath[0] = fmpdpath[0] = '\0' ;
	int    y, i, j, jmax, count;
	int    m, l, NLP, NLP0, Ntotal, N, Nclad, Nbeta, mater, profile, Nwkb, launch, Nf;
	int    Nxy;
	double lamda, k, omega, A, AA, g, n0, n1, dr, L, Tv;
	double r0, w0, dx, dy;
	double delta, NA, aa, v, w, D;
	double tau, beta, dbeta, bb, eps1, eps2, sum, Rinf, de, df, eig;
	double fmin, fmax, dfrq;					//	int    nr;
	double ncav, noxi, Avin;
	double *GI, *pulse, *q, *qg, *R, *R2, *Rb, *a, *b, *ML, *MD, *Mtau, *M, *Mbeta, *MPD;
	double **Rlp;
	int    *modem, *model, *modelin;
	int    Nvin, couple;
	int    OFFres, OFFrange;
	char   trash[65536] = "\0";
	char   directory[128] = "FILE_for_IO";
	mkdir(directory);
	/** 1. ���̓t�@�C���̓ǂݍ��� **/
	if ((fp = fopen("[FEM_Input_vcsel_guided-mode_calc].csv", "r")) != NULL) {
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

	sprintf(fp4path, "%s/FEM_result_vcsel.csv", directory);
	if ((fp4 = fopen(fp4path, "w")) != NULL) {
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
	//free_drealvector(P, 0);
	free_drealvector(M, 0);
	//
	free_drealvector(Mbeta, 0);
	free_drealvector(MPD, 0);
	free_drealvector(pulse, 0);
	//free_drealvector(opulse, 0);

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

	printf("%s\n", directory);

	//return directory[128];
}


//! //////////     �ȉ��T�u���[�`��      //////////
//
/*A.1. ��2��ό`�x�b�Z���֐� bessk (n, x) */
// ��1��ό`Bessel�֐��in=0�jI0(x)
double bessi0 (double x)
{
	double ax, ans;
	double y;
	// Polynomial fit
	if ( ( ax = fabs(x) ) < 3.75 ) {
		y = x / 3.75;
		y*= y;
		ans = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492
			+ y*(0.2659732 + y*(0.360768e-1+y*0.45813e-2)))));}
	else {
		y = 3.75 / ax;
		ans=(exp(ax) / sqrt(ax))*(0.39894228 + y*(0.1328592e-1
			+y*(0.225319e-2 + y*(-0.157565e-2 + y*(0.916281e-2
			+y*(-0.2057706e-1 + y*(0.2635537e-1 + y*(-0.1647633e-1
			+y*0.392377e-2))))))));}
	return ans;
}
// ��2��ό`Bessel�֐��in=0�j K0(x)
double bessk0 (double x)
{
	double bessi0 (double x);
	double y, ans;
	//polynomial fit
	if ( x <= 2 ) {
		y = x*x / 4.0;
		ans = (-log(x/2.0)*bessi0(x)) + (-0.57721566 + y*(0.42278420
			+ y*(0.23069756 + y*(0.3488590e-1 + y*(0.262698e-2
			+ y*(0.10750e-3 + y*0.74e-5))))));
	}
	else {
		y=2.0/x;
		ans = (exp(-x)/sqrt(x))*(1.25331414 + y*(-0.7832358e-1
			+ y*(0.2189568e-1 + y*(-0.1062446e-1 + y*(0.587872e-2
			+ y*(-0.251540e-2 + y*0.53208e-3))))));
	}
	return ans;
}
// ��1��ό`Bessel�֐��in=1�j I1(x)
double bessi1 (double x)
{
	double ax, ans;
	double y;
	if ( ( ax = fabs(x) ) < 3.75 ) {
		y = x / 3.75;
		y*=y;
		ans = ax*(0.5 + y*(0.87890594 + y*(0.51498869 + y*(0.15084934
			+ y*(0.2658733e-1 + y*(0.301532e-2 + y*0.32411e-3))))));}
	else {
		y = 3.75 / ax;
		ans = 0.2282967e-1 + y*(-0.2895312e-1 + y*(0.1787654e-1
			- y*0.420059e-2));
		ans = 0.39894228 + y*(-0.3988024e-1 + y*(-0.362018e-2
			+ y*(0.163801e-2 + y*(-0.1031555e-1 + y*ans))));
		ans*=(exp(ax) / sqrt(ax));
	}
	return x < 0.0 ? - ans : ans;
}
// ��2��ό`Bessel�֐��in=1�j K1 (x)
double bessk1 (double x) {
	double bessi1 (double x);
	double y, ans;
	if ( x <= 2.0 ) {
		y = x*x/4.0;
		ans = (log(x/2.0)*bessi1(x)) + (1.0/x)*(1.0 + y*(0.15443144
			+ y*(-0.67278579 + y*(-0.18156897 + y*(-0.1919402e-1
			+ y*(-0.110404e-2 + y*(-0.4686e-4)))))));}
	else{
		y = 2.0/x;
		ans = (exp(-x)/sqrt(x))*(1.25331414 + y*(0.23498619
			+ y*(-0.3655620e-1 + y*(0.1504268e-1 + y*(-0.780353e-2
			+ y*(0.325614e-2 + y*(-0.68245e-3)))))));
	}
	return ans;
}
// ��2��ό`Bessel�֐� Kn(x)
double bessk (int n,double x)
{
	double bessnorm = 1.0e7 ;
	double bessk0 (double x);
	double bessk1 (double x);
	int j;
	double bk,bkm,bkp,tox;
	if ( n == 0 ) return bessk0 (x) / bessnorm;
	if ( n == 1 ) return bessk1 (x) / bessnorm;
	if ( n >= 2 ) {
		tox = 2.0/x;
		bkm = bessk0(x) / bessnorm;
		bk = bessk1(x) / bessnorm;
		for ( j=1; j<n; j++ ) {
			bkp = bkm + j*tox*bk;
			bkm = bk;
			bk = bkp;
		}
		return bk;
	}
	else return 0;
}

/*A.2. ���ܗ��g�������֐�dndl_var (lamda, n_lamda, mater) */
double dndl ( double lamda, double n_lamda, int mater )
//{ return ( ( 0.01925e-4*lamda - 16.31619e-4 )*n_lamda + ( -0.02743e-4*lamda + 23.16674e-4 ) )*1.0e9; }
{ return ( ( 0.02173e-4*lamda - 18.79107e-4 )*n_lamda + ( -0.03109e-4*lamda + 26.85035e-4 ) )*1.0e9; } // 640 ~ 690 nm

/*A.3. ���ܗ��Z�x�����֐�dndl_var (lamda, mater) */ 
double dndc ( double lamda, int mater )
{ return 2.03716e-9*lamda*lamda - 3.27125e-6*lamda + 2.98314e-3; } // 589 ~ 690 nm

/*A.3. �R�A���S���ܗ� ncore (lamda, mater) */
double ncore ( double lamda, int mater ) {
	double sell; sell = 1.0;
	/* DPS-doped PMMA ( 7.8 wt.%, 1.506@589nm ) */	
	sell = sqrt ( 1.0 + ( 0.41241 / (1.0 - ( 22500 / (lamda*lamda) )) )
		+ ( 0.81215 / (1.0 - ( 6400 / (lamda*lamda) )) )
		+ ( 0.01117 / (1.0 - ( 11560000 / (lamda*lamda) )) ) );
	/* DPS-doped PMMA ( 9.0 wt.%, 1.506@655nm ) */	
//	sell = sqrt ( 1.0 + ( 0.61249 / (1.0 - ( 8467.04111 / (lamda*lamda) )) )
//		+ ( 0.61872 / (1.0 - ( 16525.45233 / (lamda*lamda) )) )
//		+ ( 0.08889 / (1.0 - ( 129601870.30589 / (lamda*lamda) )) ) );
	return sell; }

/*A.4. �N���b�h���ܗ� nclad (lamda, mater) */
double nclad ( double lamda, int mater ) {
	double sell; sell =1.0;
	/* PMMA ( 1.492@589nm ) */	
	sell = sqrt ( 1.0 + ( 0.496284 / (1.0 - ( 5154.872 / (lamda*lamda) )) )
		+ ( 0.6964977 / (1.0 - ( 13802.53 / (lamda*lamda) )) )
		+ ( 0.3223 / (1.0 - ( 85527690 / (lamda*lamda) )) ) );
	/* DPS-doped PMMA ( 1.079427062 wt. % ) */	
//	sell = sqrt ( 1.0 + ( 0.5954 / (1.0 - ( 6200.94518 / (lamda*lamda) )) )
//		+ ( 0.602 / (1.0 - ( 14730.91236 / (lamda*lamda) )) )
//		+ ( 0.29126 / (1.0 - ( 85685787.87303 / (lamda*lamda) )) ) );
	return sell; }

/*A.5. �W���s��v�f�v�Z�֐� S_matrix (a, b, q, m, n, v, w, D)*/
// b[0]=S00, b[1]=S11, ...................... b[j]=Sjj ......  b[n]=Snn
// a[0]=___, a[1]=S01, a[2]=S12, ... a[j]=Sj-1,j ... a[n]=Sn-1,n
// a[0]=___, a[1]=S10, a[2]=S21, ... a[j]=Sj,j-1 ... a[n]=Sn,n-1
void S_matrix ( double *a, double *b, double *q, int m, int n, double v, double w, double D )
{
	int j;
	if ( m == 0 ) {
			b[0] = - (1.0/2.0) + (3.0*q[0]+2.0*q[1])*(v*v/60.0)*((D*D)/(n*n)) - (1.0/12.0)*(w*w)*((D*D)/(n*n));
			b[n] = - ((2.0*n-1.0)/2.0) + ((5.0*n-2.0)*q[n-1]+3.0*(5.0*n-1.0)*q[n])*(v*v/60.0)*((D*D)/(n*n)) - ((4.0*n-1.0)/12.0)*(w*w)*((D*D)/(n*n))
				       - w*D*bessk(1,w*D)/bessk(0,w*D);}
	else {
			b[0] = 0.0;
			b[n] = - (1.0-m*m)*((2.0*n-1.0)/2.0) + ((5.0*n-2.0)*q[n-1]+3.0*(5.0*n-1.0)*q[n])*(v*v/60.0)*((D*D)/(n*n)) - ((4.0*n-1.0)/12.0)*(w*w)*((D*D)/(n*n))
				       - (m*m)*((n-1.0)*(n-1.0))*log((double)n/((double)n-1.0)) - m*m - w*D*bessk(m-1,w*D)/bessk(m,w*D) - m ; }

			b[1] = - 2.0*(1.0-m*m) + (3.0*q[0]+30*q[1]+7.0*q[2])*(v*v/60.0)*((D*D)/(n*n)) - (2.0/3.0)*(w*w)*((D*D)/(n*n)) - (m*m)*4.0*log(2.0);
		for ( j = 2; j < n; j++ )
		{	b[j] = - 2.0*j*(1.0-m*m) + ((5.0*j-2.0)*q[j-1]+30.0*j*q[j]+(5.0*j+2.0)*q[j+1])*(v*v/60.0)*((D*D)/(n*n)) - (2.0/3.0)*j*(w*w)*((D*D)/(n*n))
		              - (m*m)*( ((j-1.0)*(j-1.0))*log((double)j/((double)j-1.0)) + ((j+1.0)*(j+1.0))*log( ((double)j+1.0) / (double)j) ); }

			a[0] = 0.0;//���g�p�v�f�ɂ�0���i�[	
			a[1] = (1.0/2.0) + (2*q[0]+3*q[1])*(v*v/60.0)*((D*D)/(n*n)) - (1.0/12.0)*(w*w)*((D*D)/(n*n)) - (m*m)/2;
		for ( j = 2; j <= n; j++ )
		{	a[j] = ((2.0*j-1.0)/2.0)*(1.0-m*m) + ((5.0*j-3.0)*q[j-1] + (5.0*j-2.0)*q[j])*(v*v/60.0)*((D*D)/(n*n)) - ((2.0*j-1.0)/12.0)*(w*w)*((D*D)/(n*n))
		+ (m*m)*(j-1.0)*j*log((double)j / ((double)j-1.0)); }	
}

/*A.4. �����R���X�L�[����	mcholesky ( a, b, ML, MD, m, n )*/
void mcholesky ( double *a, double *b, double *ML, double *MD, int m, int n )
{
	int i;
	if ( m == 0 ) {
		ML[0] = 0.0;//���g�p�v�f�ɂ�0���i�[
		MD[0] = b[0];
		for ( i = 1; i <= n; i++ )
		{	MD[i] = b[i] - a[i]*a[i] / MD[i-1];
			ML[i] = a[i] / MD[i-1];	}
	}
	else {
		ML[0] = 0.0;//���g�p�v�f�ɂ�0���i�[
		ML[1] = 0.0;//���U���邩��ʈ����D0�ŗǂ��̂�?
		MD[0] = 0.0;
		MD[1] = b[1];
		for ( i = 2; i <= n; i++ )
		{	MD[i] = b[i] - a[i]*a[i] / MD[i-1];
			ML[i] = a[i] / MD[i-1];	}
	}
}

/*A.5. �����R���X�L�[����@�ɂ�������������	mcholesky_sol ( a, b, ML, MD, m, n )*/
void mcholesky_sol ( double *ML, double *MD, double *R, int m, int n )
{
	int i;
	//�uLy=R0�v������
	if ( m==0 ) {
		for ( i=1; i <= n; i++ ) {
			R[i] = R[i] - R[i-1]*ML[i]; }	
	//�u(D(LT))R1=y�v���u(LT)R1=(D-1)y=y'�v�ɕς���
		for ( i=1; i <= n; i++ ) {
			R[i] = R[i] / MD[i]; }
	//�u (LT)R1=y'�v������
		for ( i=n-1; i >= 0; i-- ) {
			R[i] = R[i] - ML[i+1]*R[i+1]; }
	}
	else{
		for ( i=2; i <= n; i++ ) {
			R[i] = R[i] - R[i-1]*ML[i]; }	
	//�u(D(LT))R1=y�v���u(LT)R1=(D-1)y=y'�v�ɕς���
		for ( i=2; i <= n; i++ ) {
			R[i] = R[i] / MD[i]; }
	//�u (LT)R1=y'�v������
		for ( i=n-1; i >= 1; i-- ) {
			R[i] = R[i] - ML[i+1]*R[i+1]; }
	}
	
}

/*A.6. �t�ׂ���@�̏����x�N�g���v�Z	R0 ( MD, R, m, n )*/
void R0 ( double *MD, double *R, int m, int n )
/* �Ίp�s��D�̐������ő�ƂȂ�v�f����1�ł���悤�ȃx�N�g����I��*/
{
	int i, j = 1;
	if ( m == 0 ) {
		for ( i = 0; i <= n-1; i++ ) {
			if ( fabs(MD[i+1]) < fabs(MD[j]) ) { j = i + 1; } }
		}
	else {
		for ( i = 1; i <= n-1; i++ ) {
			if ( fabs(MD[i+1]) < fabs(MD[j]) ) { j = i + 1; } }
	}
	//R0�̏����l�̑��
	for ( i = 0; i <= n; i++ ) {
			if ( i == j ) { R[i] = 1.0; }
			else { R[i] = 0.0; }
		}
}

/*A.7. �t�ׂ���@�̉��x�N�g���K�i��	R_norm ( R, n )*/
void R_norm ( double *R, int n )
{
	int i;
	double s = 0;
	//	�s��v�f�̂Q��a
	for ( i = 0; i <= n; i++ )	{	s = s + R[i]*R[i];	}
	if ( s	!= 0 )
	{	for ( i = 0; i <= n; i++ )	{	R[i] = R[i] / sqrt(s);	}	}
}

/*A.8. �ŗL�l�v�Z�iRayleigh quotient�jEigen ( R, a, b, m, n )*/
double Eigen ( double *R, double *a, double *b, int n, int m )
{
	// R�x�N�g�����K�i�����Ă��邽�ߓ��ς�1
	int i;
	double s=0;
	if ( m == 0 ) {
		s = ( R[0]*b[0] + R[1]*a[1] )*R[0];
		for ( i = 1; i < n; i++ ) {
			s += ( R[i-1]*a[i] + R[i]*b[i] + R[i+1]*a[i+1] )*R[i];	}
		s += ( R[n-1]*a[n] + R[n]*b[n] )*R[n];
		return s;
	}
	else {
		s = ( R[1]*b[1] + R[2]*a[2] )*R[1];
		for ( i = 2; i < n; i++ ) {
			s += ( R[i-1]*a[i] + R[i]*b[i] + R[i+1]*a[i+1] )*R[i];	}
		s += ( R[n-1]*a[n] + R[n]*b[n] )*R[n];
		return s;
	}
/*
	s = b[0]*R[0]*R[0];
	for ( i = 1; i <= n; i++ ) {
		s = s + (  b[i]*R[i]*R[i] + 2.0*a[i]*R[i-1]*R[i] );	}
*/
}

/*A.9. �Q�x���v�Z�p�֐�	dbdk_bunbo ( R, D, w, m, n )*/
/* ���̓p�����[�^�i�������d�ꕪ�z�C�R�A�a�C�`���萔�C�v�f�������j�ɑ΂���m�����[�h�Q�x���v�Z���̕��ꕪ�q���v�Z���� */

// �������d�ꐬ�� R[0]~R[n], �K�i���`���萔 w, �K�i���R�A�a D, ���ʊp���[�h���� m, ������ n
double dbdk_bunbo ( double *R, double D, double w, int m, int n )
{
	int i;
	double s=0;
	for ( i = 0; i <= n-1; i++ )
	{ s = s + (1.0/12.0) * ((D*D)/(n*n)) * ( (double)(4*i+1)*R[i]*R[i] + 2.0*(double)(2*i+1)*R[i]*R[i+1] + (double)(4*i+3)*R[i+1]*R[i+1] ); }
	if ( m == 0) {
		return s + (( bessk(1,w*D)*bessk(1,w*D) / (bessk(0,w*D)*bessk(0,w*D))) - 1.0) * ((D*D)*(R[n]*R[n]) / 2.0); }
	else {
		return s + (( bessk(m-1,w*D)*bessk(m+1,w*D) / (bessk(m,w*D)*bessk(m,w*D))) - 1.0) * ((D*D)*(R[n]*R[n]) / 2.0); }
}

/*A.10. �Q�x���v�Z�p�֐�	dbdk_bunshi ( R,q2, D, w, m, n )*/
/* ���̓p�����[�^�i�������d�ꕪ�z�C�R�A�a�C�`���萔�C�v�f�������j�ɑ΂���m�����[�h�Q�x���v�Z���̕��ꕪ�q���v�Z���� */

// �������d�ꐬ�� R[0]~R[n], �K�i���`���萔 w, �K�i���R�A�a D, ���ʊp���[�h���� m, ������ n
// ���ܗ����U�p�����[�^ q2[0]~q2[n] ( = n*(d(kn)/dk) )
double dbdk_bunshi ( double *R, double *q2, double D, double w, int m, int n)
{
	int i;
	double s=0;
	for ( i=0; i<=n-1; i++ )
	{ s = s + (1.0/12.0) * ((D*D)/(n*n))* ( ((double)(3*i)+3.0/5.0)*q2[i]*R[i]*R[i] + ((double)i+2.0/5.0)*(2.0*q2[i]*R[i+1]+q2[i+1]*R[i])*R[i] + ((double)i+3.0/5.0)*(q2[i]*R[i+1]+2.0*q2[i+1]*R[i])*R[i+1] + ((double)(3*i)+12.0/5.0)*q2[i+1]*R[i+1]*R[i+1]); }
	if ( m == 0) {
		return s + q2[n] * (( bessk(1,w*D)*bessk(1,w*D) / (bessk(0,w*D)*bessk(0,w*D))) - 1.0) * ((D*D)*(R[n]*R[n]) / 2.0); }
	else {
		return s + q2[n] * (( bessk(m-1,w*D)*bessk(m+1,w*D) / (bessk(m,w*D)*bessk(m,w*D))) - 1.0) * ((D*D)*(R[n]*R[n]) / 2.0); }
}

/* A.13. �����x�N�g���̈�m�ۗp�֐� dvector (i, j) */
int *dintvector ( int i, int j ) {
	int *a;
	if ( ( a = (int *) malloc ( (j -i+1)*sizeof (int) ) ) == NULL )
		{ printf ("Memory cannot be allocated !\n"); exit (1); }
	return (a-i); }

/* A.14. �����x�N�g���̈����p�֐� free_dvector (a, i) */
void free_dintvector ( int *a,  int i ) {
	free ( (void*) (a+i) ); }

/* A.15. �����x�N�g���̈�m�ۗp�֐� dvector (i, j) */
double *drealvector ( int i, int j ) {
	double *a;
	if ( ( a = (double *) malloc ( (j -i+1)*sizeof (double) ) ) == NULL )
		{ printf ("Memory cannot be allocated !\n"); exit (1); }
	return (a-i); }

/* A.16. �����x�N�g���̈����p�֐� free_dvector (a, i) */
void free_drealvector ( double *a,  int i ) {
	free ( (void*) (a+i) ); }

/* A.17. �����s��̈�m�ۗp�֐� dmatrix ( nr1, nr2, nl1, nl2 ) */
double **dmatrix ( int nr1, int nr2, int nl1, int nl2 ) {
	// nrow: �s�̐�, ncol: ��̐�
	double **a;
	int i, nrow, ncol;
	nrow = nr2 - nr1 +1;
	ncol  = nl2 - nl1 +1;
	/* �s�̊m�� */
	if ( ( a = (double **) malloc ( nrow*sizeof (double*) ) ) == NULL )
		{ printf ("Memory cannot be allocated !\n"); exit (1); }
	a = a - nr1; // �s�����炷
	/* ��̊m�� */
	for ( i = nr1; i <= nr2; i++ ) a[i] = (double *) malloc (ncol*sizeof (double) );
	for ( i = nr1; i <= nr2; i++ ) a[i] = a[i] - nl1;	// ������炷
	return (a); }

/* A.18.�����s��̈����p�֐� free_dmatrix ( a, nr1, nr2, nl1, nl2 ) */
void free_dmatrix ( double **a, int nr1, int nr2, int nl1, int nl2 ) {
	int i;
	for ( i = nr1; i <= nr2; i++ ) free ( (void*) (a[i] + nl1) );
	free ( (void*) (a+nr1) );
}

/* A.19. �����x�N�g���������֐� init_vector ( a, nr1, nr2 ) */
void init_intvector ( int *a, int nr1, int nr2 ) {
	int i;
	for ( i = nr1; i <= nr2; i++ ) {	a[i] = 0; }
}

/* A.20. �����x�N�g���������֐� init_vector ( a, nr1, nr2 ) */
void init_realvector ( double *a, int nr1, int nr2 ) {
	int i;
	for ( i = nr1; i <= nr2; i++ ) { a[i] = 0.0; }
}

/* A.21. �����s�񏉊����֐� init_vector ( a, nr1, nr2, nl1, nl2 ) */
void init_dmatrix (double **a, int nr1, int nr2, int nl1, int nl2 ) {
	for ( int i = nr1; i <= nr2; i++ ) {
		for ( int j = nl1; j <= nl2; j++ ) { a[i][j] = 0.0; } }
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
/* A.26. �����s�񏉊����֐� init_d3tensor ( a, nr1, nr2, nl1, nl2, np1, np2 ) */
void init_d3tensor(double ***a, int nr1, int nr2, int nl1, int nl2, int np1, int np2) {
	for (int i1 = nr1; i1 <= nr2; i1++) {
		for (int i2 = nl1; i2 <= nl2; i2++) { 
			for(int i3 = np1; i3 <= np2; i3++) { a[i1][i2][i3] = 0.0; }
}}}
/* A.24. float�^2�����z��̉�� f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh) */
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
	/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
	free((FREE_ARG)(t[nrl] + ncl - NR_END));
	free((FREE_ARG)(t + nrl - NR_END));
}
/* ��ݍ��ݐϕ� convolution (n1, P1, n2, P2) */
double* convolution ( int n1, double* P1, int n2, double* P2 ) {
	int i, j;
	double* R;
//	R = (double*) malloc ( sizeof (double) *(n1+n2+1));
	if ( ( R = (double *) malloc ( (n1+n2+1)*sizeof (double) ) ) == NULL )	
	{ printf ("Memory cannot be allocated !\n"); exit (1); }
	for ( i=0; i<=n1+n2; i++ ) { R[i] = 0.0; }
	for ( i=0; i<n1; i++ ) { 
		for( j=0; j<=n2; j++ ) { R[i+j] = R[i+j] + P1[i]*P2[j]; }}			
	//for( j=0; j<=n2; j++ ) { R[i+j]+=P1[i]*P2[j]; }}
	return R; }

/* �f�B���N�g���쐬�֐� */
void mkdir(char dirname[]) {
	struct stat statBuf;
	if (stat(dirname, &statBuf) != 0){
		if (_mkdir(dirname) != 0) {
			printf("�f�B���N�g�� %s �̍쐬�Ɏ��s���܂����B\n", dirname);
			system("pause");	}}}

//! �u�����f�B���N�g���폜��A�f�B���N�g���쐬�֐��v���s�B
/*  �u�����f�B���N�g���폜��A�f�B���N�g���쐬�֐� ���s�@�v�쐬��A�u�V���s�A�v���ȉ��̂悤�ɍ쐬�������A
	 �f�B���N�g�����t�@�C���폜�̎��s�����ǉ����ł����f�O
void delmkdirconfirm(char dirname[]) {
	struct stat statBuf;
	if (stat(dirname, &statBuf) == 0) {
		int i;
		char* str;
		sprintf(str, "rd /s %s", dirname);	//! sprintf(str, "rd /s /q %s", dirname);
		system(str);

		WIN32_FIND_DATA findData;
		HANDLE hFind;
	}
	if (_mkdir(dirname) != 0) {
		printf("�f�B���N�g�� %s �̍쐬�Ɏ��s���܂����B\n", dirname);
		system("pause");	}
}*/

//! �u�����f�B���N�g���폜��A�f�B���N�g���쐬�֐��v���s�A
/*  �u�����f�B���N�g���폜��A�f�B���N�g���쐬�֐� ���s�@�v�쐬��A�u�V���s�A�v���ȉ��̂悤�ɍ쐬�������A
     �f�B���N�g�����t�@�C���폜�̎��s�����ǉ����ł����f�O 
void delmkdir(char dirname[]) {
	struct stat statBuf;
	if (stat(dirname, &statBuf) == 0) {
		int count=0;		char listfilename[128][64];		
		TCHAR buf[256];
		char str[256];
		sprintf(str, "%s\*.csv", dirname);
		printf("%s\n", str);
#ifdef UNICODE
		MultiByteToWideChar(CP_OEMCP, MB_PRECOMPOSED, str, strlen(str), buf, (sizeof buf) / 2);
		//��������mbstowcs(buf,str,(sizeof buf)/2);
		printf("%s", buf);
#else
		strcpy(buf, str);
#endif
		WIN32_FIND_DATA findData;
		HANDLE hFind;	
		hFind = FindFirstFile(buf, &findData);
		if (hFind != INVALID_HANDLE_VALUE) {
			do {
				// �����Ƀt�@�C�����Ƃ̏������L��
				count += 1;
				printf("%d\n", count);
				//_tprintf("%s\n", findData.cFileName);	//_tprintf��printf�͓���
				sprintf(listfilename[count], "%s\%s", dirname, findData.cFileName);
				
				printf("%s\n", listfilename[count]);
			} while (FindNextFile(hFind, &findData));
			FindClose(hFind);	}
		for (int cnt = 1; cnt <= count; cnt++) {
				_rmdir(listfilename[cnt]);	}		
		_rmdir(dirname);	
	}
	if (_mkdir(dirname) != 0) {
		printf("�f�B���N�g�� %s �̍쐬�Ɏ��s���܂����B\n", dirname);
		system("pause");	}
}*/

//! �u�����f�B���N�g���폜��A�f�B���N�g���쐬�֐� ���s�@�v
//  �u�����f�B���N�g���폜��A�f�B���N�g���쐬�֐��v ���ȉ��̂悤�ɍ쐬�������A�f�B���N�g�����̃t�@�C���폜�Ɏ��s���邽�ߕύX���K�v
/*void delmkdir(char dirname[]) {
	struct stat statBuf;
	if (stat(dirname, &statBuf) == 0) {
		_rmdir(dirname);	}
	if (_mkdir(dirname) != 0) {
		printf("�f�B���N�g�� %s �̍쐬�Ɏ��s���܂����B\n", dirname);
		system("pause");	}}*/