/*
    Copyright (C) 2015, Kazutaka DOI.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "select.h"

double inner_select(int m, int n, SEXP hazard){
    if(m == 0){
        return(1);
    } else if(m > n){
        return(0);
    } else {
        return(inner_select(m, n-1, hazard) + REAL(hazard)[n-1]*inner_select(m-1, n-1, hazard));
    }
}

SEXP Rf_select(SEXP m, SEXP n, SEXP hazard){
    SEXP ret;
    m = PROTECT(coerceVector(m, INTSXP));
    n = PROTECT(coerceVector(n, INTSXP));
    hazard = PROTECT(coerceVector(hazard, REALSXP));
    ret = PROTECT(allocVector(REALSXP, 1));
    REAL(ret)[0] = inner_select(INTEGER(m)[0], INTEGER(n)[0], hazard);
    UNPROTECT(4);
    return(ret);
}
