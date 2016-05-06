/*
    Copyright (C) 2016, Kazutaka DOI.

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

#include "mumul.h"

SEXP Rf_mumul(SEXP vec){
    SEXP ret;
    vec = PROTECT(coerceVector(vec, REALSXP));
    ret = PROTECT(allocVector(REALSXP, 1));
    double val = 1.0;
    for(int i=0; i < length(vec); ++i){
        val *= REAL(vec)[i];
    }
    REAL(ret)[0] = val;
    return(ret);
}
