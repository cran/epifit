#include <cstdlib>
#include <Rinternals.h>

extern "C" SEXP Rf_scanVarname(SEXP string);

#define INT int

// Scans string vector and obtain variable name and numeric vector
// e.g., c("dose45", "dose45", "dose46") -> list("dose", c(45,46,47))
// in case of error, return list("error message", NA)
SEXP Rf_scanVarname(SEXP string) {
	int n_protect = 0;
	INT num_length = 0, charpos = 0, char_end = 0, denom_demical = 0;
	double num_begin = 0;
	int flag_numeric = 0; // one numeric value flag
	char *varname = NULL;
	const char *str;

	SEXP r_ret, r_num, r_varname;

	string = PROTECT(coerceVector(string, STRSXP));
	++n_protect;
	num_length = Rf_length(string);

	r_ret = PROTECT(allocVector(VECSXP, 2));
	++n_protect;

	r_varname = PROTECT(allocVector(STRSXP, 1));
	++n_protect;

	if (num_length < 1) {
		r_num = PROTECT(allocVector(REALSXP, 1));
		++n_protect;
		SET_STRING_ELT(r_varname, 0, Rf_mkChar("Length of argument must be larger than 0"));
		SET_VECTOR_ELT(r_ret, 0, r_varname);
		REAL(r_num)[0] = NA_REAL;
		SET_VECTOR_ELT(r_ret, 1, r_num);
		UNPROTECT(n_protect);
		return r_ret;
	}

	r_num = PROTECT(allocVector(REALSXP, num_length));
	++n_protect;

	for (INT i = 0; i < num_length; ++i) {
		charpos = 0;
		char_end = 0;
		denom_demical = 0;
		num_begin = 0;
		flag_numeric = 0;

		str = CHAR(STRING_ELT(string, i));

		if (str[0] >= '0' && str[0] <= '9') {
			r_num = PROTECT(allocVector(REALSXP, 1));
			++n_protect;
			SET_STRING_ELT(r_varname, 0, Rf_mkChar("variable name must begin character"));
			SET_VECTOR_ELT(r_ret, 0, r_varname);
			REAL(r_num)[0] = NA_REAL;
			SET_VECTOR_ELT(r_ret, 1, r_num);
			UNPROTECT(n_protect);
			return r_ret;
		}

		while (str[charpos] != '\0') {
			if (str[charpos] >= '0' && str[charpos] <= '9') {
				if (flag_numeric == 0) {
					char_end = charpos - 1;
					if (num_begin != 0) {
						r_num = PROTECT(allocVector(REALSXP, 1));
						++n_protect;
						SET_STRING_ELT(r_varname, 0, Rf_mkChar("variable name must only one numeric part"));
						SET_VECTOR_ELT(r_ret, 0, r_varname);
						REAL(r_num)[0] = NA_REAL;
						SET_VECTOR_ELT(r_ret, 1, r_num);
						UNPROTECT(n_protect);
						return r_ret;
					}
				}
				flag_numeric = 1;
				num_begin = (str[charpos] - '0') + 10 * num_begin;
				denom_demical *= 10;
			} else {
				if (flag_numeric == 1) {

					if (str[charpos] != '.') {

						r_num = PROTECT(allocVector(REALSXP, 1));
						++n_protect;
						SET_STRING_ELT(r_varname, 0, Rf_mkChar("variable end must be numeric part"));
						SET_VECTOR_ELT(r_ret, 0, r_varname);
						REAL(r_num)[0] = NA_REAL;
						SET_VECTOR_ELT(r_ret, 1, r_num);
						UNPROTECT(n_protect);
						return r_ret;

					} else {
						denom_demical = 1;
					}
				}
			}
			++charpos;
		}

		if (denom_demical != 0)
			num_begin /= denom_demical;

		if (i == 0) {

			varname = (char*)std::malloc(sizeof(char)*(char_end + 2));

			for (INT j = 0; j <= char_end; ++j)
				varname[j] = str[j];

			varname[char_end + 1] = '\0';

		} else {

			for (INT j = 0; j <= char_end; ++j)

				if (varname[j] != str[j]) {
					r_num = PROTECT(allocVector(REALSXP, 1));
					++n_protect;
					SET_STRING_ELT(r_varname, 0, Rf_mkChar("variable name is different"));
					SET_VECTOR_ELT(r_ret, 0, r_varname);
					REAL(r_num)[0] = NA_REAL;
					SET_VECTOR_ELT(r_ret, 1, r_num);
					UNPROTECT(n_protect);
					return r_ret;
				}
		}

		REAL(r_num)[i] = num_begin;
	}

	SET_STRING_ELT(r_varname, 0, Rf_mkChar(varname));
	SET_VECTOR_ELT(r_ret, 0, r_varname);
	SET_VECTOR_ELT(r_ret, 1, r_num);
	std::free(varname);
	UNPROTECT(n_protect);
	return r_ret;
}
