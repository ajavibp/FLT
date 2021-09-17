#include <flt/Kalman.hpp>

using namespace FLT;
using namespace TNT;

TNT::Array1D<double> FLT::KalmanAntec(System &Model, TNT::Array1D<double> &input, TNT::Array1D<double> &output, TNT::Array2D<double> &covariance, TNT::Array2D<double> &P, TNT::Array2D<double> &Phi)
{
	size_t n_params = Model.NumberOfAntecedents();
	size_t m = output.dim();

	if (input.dim() != Model.inputs() | m != Model.outputs())
	{
		Array1D<double> emptyV(0);
		ERRORMSGVAL(E_InNoCoherent, emptyV)
	}

	if (P.dim1() != n_params | P.dim2() !=n_params | Phi.dim1() !=n_params | Phi.dim2() !=n_params)
	{
		Array1D<double> emptyV(0);
		ERRORMSGVAL(E_NumberArgIn, emptyV)
	}

	// Preparing data
	gsl_matrix_view gsl_Phi = gsl_matrix_view_array(Phi[0], n_params, n_params);
	gsl_matrix_view gsl_P = gsl_matrix_view_array(P[0], n_params, n_params);
	Array2D<double> PPhi_t(n_params, n_params);
	gsl_matrix_view gsl_PPhi_t = gsl_matrix_view_array(PPhi_t[0], n_params, n_params);

	// A priori estimation | P = Phi*P*Phi'
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &gsl_P.matrix, &gsl_Phi.matrix, 0.0, &gsl_PPhi_t.matrix); // P * Phi'
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &gsl_Phi.matrix, &gsl_PPhi_t.matrix, 0.0, &gsl_P.matrix); // Phi * (P*Phi')

	// Estimated output
	Array1D<double> y = Model.evaluate(&input[0]);

	// Jacobian matrix of the fuzzy model with respect to its parameters
	Array2D<double> C = derantec(Model, &input[0]);

	// Filter gain
	gsl_matrix_view gsl_C = gsl_matrix_view_array(C[0], m, n_params);

	Array2D<double> temp_mn(m, n_params);
	gsl_matrix_view gsl_temp_mn = gsl_matrix_view_array(temp_mn[0], m, n_params);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &gsl_C.matrix, &gsl_P.matrix, 0.0, &gsl_temp_mn.matrix); // C * P

	// (C*P) * C'+ covariance' (A part in solve)
	gsl_matrix_view gsl_CPCt_Covart = gsl_matrix_view_array(covariance[0], m, m);
	gsl_matrix_transpose(&gsl_CPCt_Covart.matrix);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &gsl_temp_mn.matrix, &gsl_C.matrix, 1.0, &gsl_CPCt_Covart.matrix);

	// C * P*Phi' (B part in solve)
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &gsl_C.matrix, &gsl_PPhi_t.matrix, 0.0, &gsl_temp_mn.matrix);
	gsl_matrix *gsl_PhiPCt = gsl_matrix_alloc(n_params, m);
	gsl_matrix_transpose_memcpy	(gsl_PhiPCt, &gsl_temp_mn.matrix);

/*
	// ------------ Solve K' using LU (faster, but less accurate) ------------
	gsl_matrix *gsl_lu = gsl_matrix_alloc(m, m);
	gsl_matrix_memcpy(gsl_lu, &gsl_CPCt_Covart.matrix); // Solve destroys LU (for solve refine)
	gsl_permutation *perm = gsl_permutation_alloc (m);
	int signum;
	gsl_linalg_LU_decomp(gsl_lu, perm, &signum);

	Array2D<double> Kt(m,n_params,0.0);
	gsl_matrix_view gsl_Kt = gsl_matrix_view_array(Kt[0], m, n_params);

	gsl_vector_view Btemp;
	gsl_vector_view Ktemp;
	gsl_vector *residual = gsl_vector_alloc(m);
	for(size_t col=0; col<n_params; col++)
	{
		Btemp = gsl_matrix_column(&gsl_temp_mn.matrix, col);
		Ktemp = gsl_matrix_column(&gsl_Kt.matrix, col);
		gsl_linalg_LU_solve(gsl_lu, perm, &Btemp.vector, &Ktemp.vector);
		gsl_linalg_LU_refine(&gsl_CPCt_Covart.matrix, gsl_lu, perm, &Btemp.vector, &Ktemp.vector, residual);
	}
	gsl_vector_free(residual);
	gsl_permutation_free(perm);
	gsl_matrix_free(gsl_lu);
*/
	// ------------ Solve K' using QR (slower, but more accurate) ------------------------
	gsl_vector *tau = gsl_vector_alloc (m);
	gsl_linalg_QR_decomp(&gsl_CPCt_Covart.matrix, tau);

	Array2D<double> Kt(m,n_params, 0.0);
	gsl_matrix_view gsl_Kt = gsl_matrix_view_array(Kt[0], m, n_params);

	gsl_vector_view Btemp;
	gsl_vector_view Ktemp;
	gsl_vector *residual = gsl_vector_alloc(m);
	for(size_t col=0; col<n_params; col++)
	{
		Btemp = gsl_matrix_column(&gsl_temp_mn.matrix, col);
		Ktemp = gsl_matrix_column(&gsl_Kt.matrix, col);
		gsl_linalg_QR_solve(&gsl_CPCt_Covart.matrix, tau, &Btemp.vector, &Ktemp.vector);
	}
	gsl_vector_free(residual);
	gsl_vector_free(tau);

	// Parameters adjust
	Array1D<double> error = output.copy();
	gsl_vector_view gsl_error = gsl_vector_view_array(&error[0], m);
	gsl_vector_view gsl_estimated_output = gsl_vector_view_array(&y[0], m);
	gsl_vector_sub(&gsl_error.vector, &gsl_estimated_output.vector);

	Array1D<double> param = Model.getAntecedents();
	gsl_vector_view gsl_param = gsl_vector_view_array(&param[0], n_params);

	Array1D<double> temp2 = param.copy();
	gsl_vector_view gsl_temp2 = gsl_vector_view_array(&temp2[0], n_params);
	gsl_blas_dgemv(CblasNoTrans, 1.0, &gsl_Phi.matrix, &gsl_temp2.vector, 0.0, &gsl_param.vector); // Phi*param
	gsl_blas_dgemv(CblasTrans, 1.0, &gsl_Kt.matrix, &gsl_error.vector, 1.0, &gsl_param.vector); // param = (Phi*param) + K*error;

	// Model correction
	Model.setAntecedents(&param[0]);

	// Posteriori correction | P = Phi*P*Phi' - K*C*P*Phi'
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &gsl_P.matrix, &gsl_Phi.matrix, 0.0, &gsl_PPhi_t.matrix); // P * Phi'
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &gsl_Phi.matrix, &gsl_PPhi_t.matrix, 0.0, &gsl_P.matrix); // Phi * (P*Phi')
	gsl_blas_dgemm(CblasTrans, CblasTrans, -1.0, &gsl_Kt.matrix, gsl_PhiPCt, 1.0, &gsl_P.matrix); // P = Phi*P*Phi' - K*(C*P*PHi')

	gsl_matrix_free(gsl_PhiPCt);
	return error;
}

TNT::Array1D<double> FLT::KalmanAntec(System &Model, TNT::Array1D<double> &input, TNT::Array1D<double> &output, TNT::Array2D<double> &covariance, TNT::Array2D<double> &P)
{
	size_t n_params = Model.NumberOfAntecedents();
	size_t m = output.dim();

	if (input.dim() != Model.inputs() | m != Model.outputs())
	{
		Array1D<double> emptyV(0);
		ERRORMSGVAL(E_InNoCoherent, emptyV)
	}

	if (P.dim1() != n_params | P.dim2() !=n_params)
	{
		Array1D<double> emptyV(0);
		ERRORMSGVAL(E_NumberArgIn, emptyV)
	}

	// Preparing data
	gsl_matrix_view gsl_P = gsl_matrix_view_array(P[0], n_params, n_params);

	// A priori estimation it is not necessary, P = Phi*P*Phi' where Phi is the identity matrix

	// Estimated output
	Array1D<double> y = Model.evaluate(&input[0]);

	// Jacobian matrix of the fuzzy model with respect to its parameters
	Array2D<double> C = derantec(Model, &input[0]);

	// Filter gain
	gsl_matrix_view gsl_C = gsl_matrix_view_array(C[0], m, n_params);

	Array2D<double> temp_mn(m, n_params);
	gsl_matrix_view gsl_temp_mn = gsl_matrix_view_array(temp_mn[0], m, n_params);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &gsl_C.matrix, &gsl_P.matrix, 0.0, &gsl_temp_mn.matrix); // C * P

	// (C*P) * C'+ covariance' (A part in solve)
	gsl_matrix_view gsl_CPCt_Covart = gsl_matrix_view_array(covariance[0], m, m);
	gsl_matrix_transpose(&gsl_CPCt_Covart.matrix);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &gsl_temp_mn.matrix, &gsl_C.matrix, 1.0, &gsl_CPCt_Covart.matrix);

	// C * P*Phi' (B part in solve), Phi is the identity matrix
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &gsl_C.matrix, &gsl_P.matrix, 0.0, &gsl_temp_mn.matrix);
	gsl_matrix *gsl_PhiPCt = gsl_matrix_alloc(n_params, m);
	gsl_matrix_transpose_memcpy	(gsl_PhiPCt, &gsl_temp_mn.matrix);

/*
	// ------------ Solve K' using LU (faster but less accurate) ------------
	gsl_matrix *gsl_lu = gsl_matrix_alloc(m, m);
	gsl_matrix_memcpy(gsl_lu, &gsl_CPCt_Covart.matrix); // Solve destroys LU (for solve refine)
	gsl_permutation *perm = gsl_permutation_alloc (m);
	int signum;
	gsl_linalg_LU_decomp(gsl_lu, perm, &signum);

	Array2D<double> Kt(m,n_params,0.0);
	gsl_matrix_view gsl_Kt = gsl_matrix_view_array(Kt[0], m, n_params);

	gsl_vector_view Btemp;
	gsl_vector_view Ktemp;
	gsl_vector *residual = gsl_vector_alloc(m);
	for(size_t col=0; col<n_params; col++)
	{
		Btemp = gsl_matrix_column(&gsl_temp_mn.matrix, col);
		Ktemp = gsl_matrix_column(&gsl_Kt.matrix, col);
		gsl_linalg_LU_solve(gsl_lu, perm, &Btemp.vector, &Ktemp.vector);
		gsl_linalg_LU_refine(&gsl_CPCt_Covart.matrix, gsl_lu, perm, &Btemp.vector, &Ktemp.vector, residual);
	}
	gsl_vector_free(residual);
	gsl_permutation_free(perm);
	gsl_matrix_free(gsl_lu);
*/
	// ------------ Solve K' using QR ------------
	gsl_vector *tau = gsl_vector_alloc (m);
	gsl_linalg_QR_decomp(&gsl_CPCt_Covart.matrix, tau);

	Array2D<double> Kt(m,n_params, 0.0);
	gsl_matrix_view gsl_Kt = gsl_matrix_view_array(Kt[0], m, n_params);

	gsl_vector_view Btemp;
	gsl_vector_view Ktemp;
	gsl_vector *residual = gsl_vector_alloc(m);
	for(size_t col=0; col<n_params; col++)
	{
		Btemp = gsl_matrix_column(&gsl_temp_mn.matrix, col);
		Ktemp = gsl_matrix_column(&gsl_Kt.matrix, col);
		gsl_linalg_QR_solve(&gsl_CPCt_Covart.matrix, tau, &Btemp.vector, &Ktemp.vector);
	}
	gsl_vector_free(residual);
	gsl_vector_free(tau);

	// Parameters adjust
	Array1D<double> error = output.copy();
	gsl_vector_view gsl_error = gsl_vector_view_array(&error[0], m);
	gsl_vector_view gsl_estimated_output = gsl_vector_view_array(&y[0], m);
	gsl_vector_sub(&gsl_error.vector, &gsl_estimated_output.vector);

	Array1D<double> param = Model.getAntecedents();
	gsl_vector_view gsl_param = gsl_vector_view_array(&param[0], n_params);

	gsl_blas_dgemv(CblasTrans, 1.0, &gsl_Kt.matrix, &gsl_error.vector, 1.0, &gsl_param.vector); // param = (Phi*param) + K*error; Phi is the identity matrix

	// Model correction
	Model.setAntecedents(&param[0]);

	// Posteriori correction | P = Phi*P*Phi' - K*C*P*Phi'
	gsl_blas_dgemm(CblasTrans, CblasTrans, -1.0, &gsl_Kt.matrix, gsl_PhiPCt, 1.0, &gsl_P.matrix); // P = Phi*P*Phi' - K*(C*P*PHi')

	gsl_matrix_free(gsl_PhiPCt);
	return error;
}

TNT::Array1D<double> FLT::KalmanConseq(System &Model, TNT::Array1D<double> &input, TNT::Array1D<double> &output, TNT::Array2D<double> &covariance, TNT::Array2D<double> &P, TNT::Array2D<double> &Phi)
{
	size_t n_params = Model.NumberOfConsequents();
	size_t m = output.dim();

	if (input.dim() != Model.inputs() | m != Model.outputs())
	{
		Array1D<double> emptyV(0);
		ERRORMSGVAL(E_InNoCoherent, emptyV)
	}

	if (P.dim1() != n_params | P.dim2() !=n_params | Phi.dim1() !=n_params | Phi.dim2() !=n_params)
	{
		Array1D<double> emptyV(0);
		ERRORMSGVAL(E_NumberArgIn, emptyV)
	}

	// Preparing data
	gsl_matrix_view gsl_Phi = gsl_matrix_view_array(Phi[0], n_params, n_params);
	gsl_matrix_view gsl_P = gsl_matrix_view_array(P[0], n_params, n_params);
	Array2D<double> PPhi_t(n_params, n_params);
	gsl_matrix_view gsl_PPhi_t = gsl_matrix_view_array(PPhi_t[0], n_params, n_params);

	// A priori estimation | P = Phi*P*Phi'
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &gsl_P.matrix, &gsl_Phi.matrix, 0.0, &gsl_PPhi_t.matrix); // P * Phi'
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &gsl_Phi.matrix, &gsl_PPhi_t.matrix, 0.0, &gsl_P.matrix); // Phi * (P*Phi')

	// Estimated output
	Array1D<double> y = Model.evaluate(&input[0]);

	// Jacobian matrix of the fuzzy model with respect to its parameters
	Array2D<double> C = derconseq(Model, &input[0]);

	// Filter gain
	gsl_matrix_view gsl_C = gsl_matrix_view_array(C[0], m, n_params);

	Array2D<double> temp_mn(m, n_params);
	gsl_matrix_view gsl_temp_mn = gsl_matrix_view_array(temp_mn[0], m, n_params);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &gsl_C.matrix, &gsl_P.matrix, 0.0, &gsl_temp_mn.matrix); // C * P

	// (C*P) * C'+ covariance' (A part in solve)
	gsl_matrix_view gsl_CPCt_Covart = gsl_matrix_view_array(covariance[0], m, m);
	gsl_matrix_transpose(&gsl_CPCt_Covart.matrix);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &gsl_temp_mn.matrix, &gsl_C.matrix, 1.0, &gsl_CPCt_Covart.matrix);

	// C * P*Phi' (B part in solve)
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &gsl_C.matrix, &gsl_PPhi_t.matrix, 0.0, &gsl_temp_mn.matrix);
	gsl_matrix *gsl_PhiPCt = gsl_matrix_alloc(n_params, m);
	gsl_matrix_transpose_memcpy	(gsl_PhiPCt, &gsl_temp_mn.matrix);

/*
	// ------------ Solve K' using LU (faster but less accurate) ------------
	gsl_matrix *gsl_lu = gsl_matrix_alloc(m, m);
	gsl_matrix_memcpy(gsl_lu, &gsl_CPCt_Covart.matrix); // Solve destroys LU (for solve refine)
	gsl_permutation *perm = gsl_permutation_alloc (m);
	int signum;
	gsl_linalg_LU_decomp(gsl_lu, perm, &signum);

	Array2D<double> Kt(m,n_params,0.0);
	gsl_matrix_view gsl_Kt = gsl_matrix_view_array(Kt[0], m, n_params);

	gsl_vector_view Btemp;
	gsl_vector_view Ktemp;
	gsl_vector *residual = gsl_vector_alloc(m);
	for(size_t col=0; col<n_params; col++)
	{
		Btemp = gsl_matrix_column(&gsl_temp_mn.matrix, col);
		Ktemp = gsl_matrix_column(&gsl_Kt.matrix, col);
		gsl_linalg_LU_solve(gsl_lu, perm, &Btemp.vector, &Ktemp.vector);
		gsl_linalg_LU_refine(&gsl_CPCt_Covart.matrix, gsl_lu, perm, &Btemp.vector, &Ktemp.vector, residual);
	}
	gsl_vector_free(residual);
	gsl_permutation_free(perm);
	gsl_matrix_free(gsl_lu);
*/
	// ------------ Solve K' using QR ------------
	gsl_vector *tau = gsl_vector_alloc (m);
	gsl_linalg_QR_decomp(&gsl_CPCt_Covart.matrix, tau);

	Array2D<double> Kt(m,n_params, 0.0);
	gsl_matrix_view gsl_Kt = gsl_matrix_view_array(Kt[0], m, n_params);

	gsl_vector_view Btemp;
	gsl_vector_view Ktemp;
	gsl_vector *residual = gsl_vector_alloc(m);
	for(size_t col=0; col<n_params; col++)
	{
		Btemp = gsl_matrix_column(&gsl_temp_mn.matrix, col);
		Ktemp = gsl_matrix_column(&gsl_Kt.matrix, col);
		gsl_linalg_QR_solve(&gsl_CPCt_Covart.matrix, tau, &Btemp.vector, &Ktemp.vector);
	}
	gsl_vector_free(residual);
	gsl_vector_free(tau);

	// Parameters adjust
	Array1D<double> error = output.copy();
	gsl_vector_view gsl_error = gsl_vector_view_array(&error[0], m);
	gsl_vector_view gsl_estimated_output = gsl_vector_view_array(&y[0], m);
	gsl_vector_sub(&gsl_error.vector, &gsl_estimated_output.vector);

	Array1D<double> param = Model.getConsequents();
	gsl_vector_view gsl_param = gsl_vector_view_array(&param[0], n_params);

	Array1D<double> temp2 = param.copy();
	gsl_vector_view gsl_temp2 = gsl_vector_view_array(&temp2[0], n_params);
	gsl_blas_dgemv(CblasNoTrans, 1.0, &gsl_Phi.matrix, &gsl_temp2.vector, 0.0, &gsl_param.vector); // Phi*param
	gsl_blas_dgemv(CblasTrans, 1.0, &gsl_Kt.matrix, &gsl_error.vector, 1.0, &gsl_param.vector); // param = (Phi*param) + K*error;

	// Model correction
	Model.setConsequents(&param[0]);

	// Posteriori correction | P = Phi*P*Phi' - K*C*P*Phi'
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &gsl_P.matrix, &gsl_Phi.matrix, 0.0, &gsl_PPhi_t.matrix); // P * Phi'
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &gsl_Phi.matrix, &gsl_PPhi_t.matrix, 0.0, &gsl_P.matrix); // Phi * (P*Phi')
	gsl_blas_dgemm(CblasTrans, CblasTrans, -1.0, &gsl_Kt.matrix, gsl_PhiPCt, 1.0, &gsl_P.matrix); // P = Phi*P*Phi' - K*(C*P*PHi')

	gsl_matrix_free(gsl_PhiPCt);
	return error;
}

TNT::Array1D<double> FLT::KalmanConseq(System &Model, TNT::Array1D<double> &input, TNT::Array1D<double> &output, TNT::Array2D<double> &covariance, TNT::Array2D<double> &P)
{
	size_t n_params = Model.NumberOfConsequents();
	size_t m = output.dim();

	if (input.dim() != Model.inputs() | m != Model.outputs())
	{
		Array1D<double> emptyV(0);
		ERRORMSGVAL(E_InNoCoherent, emptyV)
	}

	if (P.dim1() != n_params | P.dim2() !=n_params)
	{
		Array1D<double> emptyV(0);
		ERRORMSGVAL(E_NumberArgIn, emptyV)
	}

	// Preparing data
	gsl_matrix_view gsl_P = gsl_matrix_view_array(P[0], n_params, n_params);

	// A priori estimation it is not necessary, P = Phi*P*Phi' where Phi is the identity matrix

	// Estimated output
	Array1D<double> y = Model.evaluate(&input[0]);

	// Jacobian matrix of the fuzzy model with respect to its parameters
	Array2D<double> C = derconseq(Model, &input[0]);

	// Filter gain
	gsl_matrix_view gsl_C = gsl_matrix_view_array(C[0], m, n_params);

	Array2D<double> temp_mn(m, n_params);
	gsl_matrix_view gsl_temp_mn = gsl_matrix_view_array(temp_mn[0], m, n_params);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &gsl_C.matrix, &gsl_P.matrix, 0.0, &gsl_temp_mn.matrix); // C * P

	// (C*P) * C'+ covariance' (A part in solve)
	gsl_matrix_view gsl_CPCt_Covart = gsl_matrix_view_array(covariance[0], m, m);
	gsl_matrix_transpose(&gsl_CPCt_Covart.matrix);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &gsl_temp_mn.matrix, &gsl_C.matrix, 1.0, &gsl_CPCt_Covart.matrix);

	// C * P*Phi' (B part in solve)
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &gsl_C.matrix, &gsl_P.matrix, 0.0, &gsl_temp_mn.matrix);
	gsl_matrix *gsl_PhiPCt = gsl_matrix_alloc(n_params, m);
	gsl_matrix_transpose_memcpy	(gsl_PhiPCt, &gsl_temp_mn.matrix);

/*
	// ------------ Solve K' using LU (faster but less accurate) ------------
	gsl_matrix *gsl_lu = gsl_matrix_alloc(m, m);
	gsl_matrix_memcpy(gsl_lu, &gsl_CPCt_Covart.matrix); // Solve destroys LU (for solve refine)
	gsl_permutation *perm = gsl_permutation_alloc (m);
	int signum;
	gsl_linalg_LU_decomp(gsl_lu, perm, &signum);

	Array2D<double> Kt(m,n_params,0.0);
	gsl_matrix_view gsl_Kt = gsl_matrix_view_array(Kt[0], m, n_params);

	gsl_vector_view Btemp;
	gsl_vector_view Ktemp;
	gsl_vector *residual = gsl_vector_alloc(m);
	for(size_t col=0; col<n_params; col++)
	{
		Btemp = gsl_matrix_column(&gsl_temp_mn.matrix, col);
		Ktemp = gsl_matrix_column(&gsl_Kt.matrix, col);
		gsl_linalg_LU_solve(gsl_lu, perm, &Btemp.vector, &Ktemp.vector);
		gsl_linalg_LU_refine(&gsl_CPCt_Covart.matrix, gsl_lu, perm, &Btemp.vector, &Ktemp.vector, residual);
	}
	gsl_vector_free(residual);
	gsl_permutation_free(perm);
	gsl_matrix_free(gsl_lu);
*/
	// ------------ Solve K' using QR ------------
	gsl_vector *tau = gsl_vector_alloc (m);
	gsl_linalg_QR_decomp(&gsl_CPCt_Covart.matrix, tau);

	Array2D<double> Kt(m,n_params, 0.0);
	gsl_matrix_view gsl_Kt = gsl_matrix_view_array(Kt[0], m, n_params);

	gsl_vector_view Btemp;
	gsl_vector_view Ktemp;
	gsl_vector *residual = gsl_vector_alloc(m);
	for(size_t col=0; col<n_params; col++)
	{
		Btemp = gsl_matrix_column(&gsl_temp_mn.matrix, col);
		Ktemp = gsl_matrix_column(&gsl_Kt.matrix, col);
		gsl_linalg_QR_solve(&gsl_CPCt_Covart.matrix, tau, &Btemp.vector, &Ktemp.vector);
	}
	gsl_vector_free(residual);
	gsl_vector_free(tau);

	// Parameters adjust
	Array1D<double> error = output.copy();
	gsl_vector_view gsl_error = gsl_vector_view_array(&error[0], m);
	gsl_vector_view gsl_estimated_output = gsl_vector_view_array(&y[0], m);
	gsl_vector_sub(&gsl_error.vector, &gsl_estimated_output.vector);

	Array1D<double> param = Model.getConsequents();
	gsl_vector_view gsl_param = gsl_vector_view_array(&param[0], n_params);

	gsl_blas_dgemv(CblasTrans, 1.0, &gsl_Kt.matrix, &gsl_error.vector, 1.0, &gsl_param.vector); // param = (Phi*param) + K*error; Phi is the identity matrix

	// Model correction
	Model.setConsequents(&param[0]);

	// Posteriori correction | P = Phi*P*Phi' - K*C*P*Phi', Phi is the identity matrix
	gsl_blas_dgemm(CblasTrans, CblasTrans, -1.0, &gsl_Kt.matrix, gsl_PhiPCt, 1.0, &gsl_P.matrix); // P = Phi*P*Phi' - K*(C*P*PHi')

	gsl_matrix_free(gsl_PhiPCt);
	return error;
}

TNT::Array1D<double> FLT::KalmanFuz(System &Model, TNT::Array1D<double> &input, TNT::Array1D<double> &output, TNT::Array2D<double> &covariance, TNT::Array2D<double> &P, TNT::Array2D<double> &Phi)
{
	size_t num_antecedents = Model.NumberOfAntecedents();
	size_t n_params = Model.NumberOfConsequents() + num_antecedents;
	size_t m = output.dim();

	if (input.dim() != Model.inputs() | m != Model.outputs())
	{
		Array1D<double> emptyV(0);
		ERRORMSGVAL(E_InNoCoherent, emptyV)
	}

	if (P.dim1() != n_params | P.dim2() !=n_params | Phi.dim1() !=n_params | Phi.dim2() !=n_params)
	{
		Array1D<double> emptyV(0);
		ERRORMSGVAL(E_NumberArgIn, emptyV)
	}

	// Preparing data
	gsl_matrix_view gsl_Phi = gsl_matrix_view_array(Phi[0], n_params, n_params);
	gsl_matrix_view gsl_P = gsl_matrix_view_array(P[0], n_params, n_params);
	Array2D<double> PPhi_t(n_params, n_params);
	gsl_matrix_view gsl_PPhi_t = gsl_matrix_view_array(PPhi_t[0], n_params, n_params);

	// A priori estimation | P = Phi*P*Phi'
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &gsl_P.matrix, &gsl_Phi.matrix, 0.0, &gsl_PPhi_t.matrix); // P * Phi'
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &gsl_Phi.matrix, &gsl_PPhi_t.matrix, 0.0, &gsl_P.matrix); // Phi * (P*Phi')

	// Estimated output
	Array1D<double> y = Model.evaluate(&input[0]);

	// Jacobian matrix of the fuzzy model with respect to its parameters
	Array2D<double> C = derfuzzy(Model, &input[0]);

	// Filter gain
	gsl_matrix_view gsl_C = gsl_matrix_view_array(C[0], m, n_params);

	Array2D<double> temp_mn(m, n_params);
	gsl_matrix_view gsl_temp_mn = gsl_matrix_view_array(temp_mn[0], m, n_params);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &gsl_C.matrix, &gsl_P.matrix, 0.0, &gsl_temp_mn.matrix); // C * P

	// (C*P) * C'+ covariance' (A part in solve)
	gsl_matrix_view gsl_CPCt_Covart = gsl_matrix_view_array(covariance[0], m, m);
	gsl_matrix_transpose(&gsl_CPCt_Covart.matrix);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &gsl_temp_mn.matrix, &gsl_C.matrix, 1.0, &gsl_CPCt_Covart.matrix);

	// C * P*Phi' (B part in solve)
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &gsl_C.matrix, &gsl_PPhi_t.matrix, 0.0, &gsl_temp_mn.matrix);
	gsl_matrix *gsl_PhiPCt = gsl_matrix_alloc(n_params, m);
	gsl_matrix_transpose_memcpy	(gsl_PhiPCt, &gsl_temp_mn.matrix);

/*
	// ------------ Solve K' using LU (faster but less accurate) ------------
	gsl_matrix *gsl_lu = gsl_matrix_alloc(m, m);
	gsl_matrix_memcpy(gsl_lu, &gsl_CPCt_Covart.matrix); // Solve destroys LU (for solve refine)
	gsl_permutation *perm = gsl_permutation_alloc (m);
	int signum;
	gsl_linalg_LU_decomp(gsl_lu, perm, &signum);

	Array2D<double> Kt(m,n_params,0.0);
	gsl_matrix_view gsl_Kt = gsl_matrix_view_array(Kt[0], m, n_params);

	gsl_vector_view Btemp;
	gsl_vector_view Ktemp;
	gsl_vector *residual = gsl_vector_alloc(m);
	for(size_t col=0; col<n_params; col++)
	{
		Btemp = gsl_matrix_column(&gsl_temp_mn.matrix, col);
		Ktemp = gsl_matrix_column(&gsl_Kt.matrix, col);
		gsl_linalg_LU_solve(gsl_lu, perm, &Btemp.vector, &Ktemp.vector);
		gsl_linalg_LU_refine(&gsl_CPCt_Covart.matrix, gsl_lu, perm, &Btemp.vector, &Ktemp.vector, residual);
	}
	gsl_vector_free(residual);
	gsl_permutation_free(perm);
	gsl_matrix_free(gsl_lu);
*/
	// ------------ Solve K' using QR ------------
	gsl_vector *tau = gsl_vector_alloc (m);
	gsl_linalg_QR_decomp(&gsl_CPCt_Covart.matrix, tau);

	Array2D<double> Kt(m,n_params, 0.0);
	gsl_matrix_view gsl_Kt = gsl_matrix_view_array(Kt[0], m, n_params);

	gsl_vector_view Btemp;
	gsl_vector_view Ktemp;
	gsl_vector *residual = gsl_vector_alloc(m);
	for(size_t col=0; col<n_params; col++)
	{
		Btemp = gsl_matrix_column(&gsl_temp_mn.matrix, col);
		Ktemp = gsl_matrix_column(&gsl_Kt.matrix, col);
		gsl_linalg_QR_solve(&gsl_CPCt_Covart.matrix, tau, &Btemp.vector, &Ktemp.vector);
	}
	gsl_vector_free(residual);
	gsl_vector_free(tau);

	// Parameters adjust
	Array1D<double> error = output.copy();
	gsl_vector_view gsl_error = gsl_vector_view_array(&error[0], m);
	gsl_vector_view gsl_estimated_output = gsl_vector_view_array(&y[0], m);
	gsl_vector_sub(&gsl_error.vector, &gsl_estimated_output.vector);

	Array1D<double> param(n_params);
	Array1D<double> AntecData = Model.getAntecedents();
	Array1D<double> ConseqData = Model.getConsequents();
	double *p_param = &param[0];
	for (size_t d=0;d<AntecData.dim();d++, p_param++)
		*p_param = AntecData[d];
	for (size_t d=0;d<ConseqData.dim();d++, p_param++)
		*p_param = ConseqData[d];

	gsl_vector_view gsl_param = gsl_vector_view_array(&param[0], n_params);

	Array1D<double> temp2 = param.copy();
	gsl_vector_view gsl_temp2 = gsl_vector_view_array(&temp2[0], n_params);
	gsl_blas_dgemv(CblasNoTrans, 1.0, &gsl_Phi.matrix, &gsl_temp2.vector, 0.0, &gsl_param.vector); // Phi*param
	gsl_blas_dgemv(CblasTrans, 1.0, &gsl_Kt.matrix, &gsl_error.vector, 1.0, &gsl_param.vector); // param = (Phi*param) + K*error;

	// Model correction
	p_param = &param[0];
	Model.setAntecedents(p_param);
	p_param += num_antecedents;
	Model.setConsequents(p_param);

	// Posteriori correction | P = Phi*P*Phi' - K*C*P*Phi'
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &gsl_P.matrix, &gsl_Phi.matrix, 0.0, &gsl_PPhi_t.matrix); // P * Phi'
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &gsl_Phi.matrix, &gsl_PPhi_t.matrix, 0.0, &gsl_P.matrix); // Phi * (P*Phi')
	gsl_blas_dgemm(CblasTrans, CblasTrans, -1.0, &gsl_Kt.matrix, gsl_PhiPCt, 1.0, &gsl_P.matrix); // P = Phi*P*Phi' - K*(C*P*PHi')

	gsl_matrix_free(gsl_PhiPCt);
	return error;
}

TNT::Array1D<double> FLT::KalmanFuz(System &Model, TNT::Array1D<double> &input, TNT::Array1D<double> &output, TNT::Array2D<double> &covariance, TNT::Array2D<double> &P)
{
	size_t num_antecedents = Model.NumberOfAntecedents();
	size_t n_params = Model.NumberOfConsequents() + num_antecedents;
	size_t m = output.dim();

	if (input.dim() != Model.inputs() | m != Model.outputs())
	{
		Array1D<double> emptyV(0);
		ERRORMSGVAL(E_InNoCoherent, emptyV)
	}

	if (P.dim1() != n_params | P.dim2() !=n_params)
	{
		Array1D<double> emptyV(0);
		ERRORMSGVAL(E_NumberArgIn, emptyV)
	}

	// Preparing data
	gsl_matrix_view gsl_P = gsl_matrix_view_array(P[0], n_params, n_params);

	// A priori estimation it is not necessary, P = Phi*P*Phi' where Phi is the identity matrix

	// Estimated output
	Array1D<double> y = Model.evaluate(&input[0]);

	// Jacobian matrix of the fuzzy model with respect to its parameters
	Array2D<double> C = derfuzzy(Model, &input[0]);

	// Filter gain
	gsl_matrix_view gsl_C = gsl_matrix_view_array(C[0], m, n_params);

	Array2D<double> temp_mn(m, n_params);
	gsl_matrix_view gsl_temp_mn = gsl_matrix_view_array(temp_mn[0], m, n_params);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &gsl_C.matrix, &gsl_P.matrix, 0.0, &gsl_temp_mn.matrix); // C * P

	// (C*P) * C'+ covariance' (A part in solve)
	gsl_matrix_view gsl_CPCt_Covart = gsl_matrix_view_array(covariance[0], m, m);
	gsl_matrix_transpose(&gsl_CPCt_Covart.matrix);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &gsl_temp_mn.matrix, &gsl_C.matrix, 1.0, &gsl_CPCt_Covart.matrix);

	// C * P*Phi' (B part in solve), Phi is the identity matrix
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &gsl_C.matrix, &gsl_P.matrix, 0.0, &gsl_temp_mn.matrix);
	gsl_matrix *gsl_PhiPCt = gsl_matrix_alloc(n_params, m);
	gsl_matrix_transpose_memcpy	(gsl_PhiPCt, &gsl_temp_mn.matrix);

/*
	// ------------ Solve K' using LU (faster but less accurate) ------------
	gsl_matrix *gsl_lu = gsl_matrix_alloc(m, m);
	gsl_matrix_memcpy(gsl_lu, &gsl_CPCt_Covart.matrix); // Solve destroys LU (for solve refine)
	gsl_permutation *perm = gsl_permutation_alloc (m);
	int signum;
	gsl_linalg_LU_decomp(gsl_lu, perm, &signum);

	Array2D<double> Kt(m,n_params,0.0);
	gsl_matrix_view gsl_Kt = gsl_matrix_view_array(Kt[0], m, n_params);

	gsl_vector_view Btemp;
	gsl_vector_view Ktemp;
	gsl_vector *residual = gsl_vector_alloc(m);
	for(size_t col=0; col<n_params; col++)
	{
		Btemp = gsl_matrix_column(&gsl_temp_mn.matrix, col);
		Ktemp = gsl_matrix_column(&gsl_Kt.matrix, col);
		gsl_linalg_LU_solve(gsl_lu, perm, &Btemp.vector, &Ktemp.vector);
		gsl_linalg_LU_refine(&gsl_CPCt_Covart.matrix, gsl_lu, perm, &Btemp.vector, &Ktemp.vector, residual);
	}
	gsl_vector_free(residual);
	gsl_permutation_free(perm);
	gsl_matrix_free(gsl_lu);
*/
	// ------------ Solve K' using QR ------------
	gsl_vector *tau = gsl_vector_alloc (m);
	gsl_linalg_QR_decomp(&gsl_CPCt_Covart.matrix, tau);

	Array2D<double> Kt(m,n_params, 0.0);
	gsl_matrix_view gsl_Kt = gsl_matrix_view_array(Kt[0], m, n_params);

	gsl_vector_view Btemp;
	gsl_vector_view Ktemp;
	gsl_vector *residual = gsl_vector_alloc(m);
	for(size_t col=0; col<n_params; col++)
	{
		Btemp = gsl_matrix_column(&gsl_temp_mn.matrix, col);
		Ktemp = gsl_matrix_column(&gsl_Kt.matrix, col);
		gsl_linalg_QR_solve(&gsl_CPCt_Covart.matrix, tau, &Btemp.vector, &Ktemp.vector);
	}
	gsl_vector_free(residual);
	gsl_vector_free(tau);

	// Parameters adjust
	Array1D<double> error = output.copy();
	gsl_vector_view gsl_error = gsl_vector_view_array(&error[0], m);
	gsl_vector_view gsl_estimated_output = gsl_vector_view_array(&y[0], m);
	gsl_vector_sub(&gsl_error.vector, &gsl_estimated_output.vector);

	Array1D<double> param(n_params);
	Array1D<double> AntecData = Model.getAntecedents();
	Array1D<double> ConseqData = Model.getConsequents();
	double *p_param = &param[0];
	for (size_t d=0;d<AntecData.dim();d++, p_param++)
		*p_param = AntecData[d];
	for (size_t d=0;d<ConseqData.dim();d++, p_param++)
		*p_param = ConseqData[d];

	gsl_vector_view gsl_param = gsl_vector_view_array(&param[0], n_params);

	gsl_blas_dgemv(CblasTrans, 1.0, &gsl_Kt.matrix, &gsl_error.vector, 1.0, &gsl_param.vector); // param = (Phi*param) + K*error; , Phi is the identity matrix

	// Model correction
	p_param = &param[0];
	Model.setAntecedents(p_param);
	p_param += num_antecedents;
	Model.setConsequents(p_param);

	// Posteriori correction | P = Phi*P*Phi' - K*C*P*Phi', Phi is the identity matrix
	gsl_blas_dgemm(CblasTrans, CblasTrans, -1.0, &gsl_Kt.matrix, gsl_PhiPCt, 1.0, &gsl_P.matrix); // P = Phi*P*Phi' - K*(C*P*PHi'), Phi is the identity matrix

	gsl_matrix_free(gsl_PhiPCt);
	return error;
}
