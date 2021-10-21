/* tanimoto: Module for pairwise tanimoto calculations in C from python
 * 
 * Function List:
 *
 * bitcounts(fplist)
 * Calculate number of bits on in raw binary fingerprints
 *
 * fractionfilterB(fraction, fplistA, countsA, fplistB, countsB)
 * Return top fraction of listB indices ordered by 
 *     decreasing maximum tanimoto to listA
 *
 * simfilterB(cutoff, fplistA, countsA, fplistB, countsB)
 * Return indices of listB that are more similar than cutoff to listA
 *
 * matrix(fplistA, countsA, fplistB, countsB)
 * Calculate matrix of tanimoto coefficients into python list of lists.
 *
 * Revisions:
 * Mysinger 7/05 Created
 * Mysinger 8/05 Modified for new dud
 * Mysinger 5/07 Modified for Py_ssize_t
 * Mysinger 9/11 Modified for fractionfilterB
 */

#include <Python.h>
#include <stdlib.h>
#include "fast_tanimoto.h"

#if PY_VERSION_HEX < 0x02050000 && !defined(PY_SSIZE_T_MIN)
typedef int Py_ssize_t;
#define PY_SSIZE_T_MAX INT_MAX
#define PY_SSIZE_T_MIN INT_MIN
#endif

#define FreeAB { free(fpA); free(fpB); free(countA); free(countB); }

static Py_ssize_t copy_fp_list(PyObject *pylist, Py_ssize_t pylen, 
			       uint16_t **clist)
{
  PyObject *x;
  char *fptemp;
  Py_ssize_t fplen, dblen, i;

  /* Move first list element to get the fingerprint double byte length */
  x = PyList_GetItem(pylist, 0);
  if (PyString_AsStringAndSize(x, &fptemp, &fplen) == -1) {
    return -1;
  }
  clist[0] = (uint16_t *) fptemp;
  dblen = fplen/2;

  /* Move list from python list to c array */
  for (i = 1; i < pylen; i++) {
    x = PyList_GetItem(pylist, i);
    if (PyString_AsStringAndSize(x, &fptemp, &fplen) == -1) {
      return -1;
    }
    clist[i] = (uint16_t *) fptemp;
    if ((dblen*2) != fplen) {
      PyErr_SetString(PyExc_ValueError, 
		      "All fingerprints must be the same length!");
      return -1;
    }
  }

  return dblen;
}

static int copy_count_list(PyObject *pylist, Py_ssize_t pylen, int *clist)
{
  Py_ssize_t i;

  /* Move list from python list to c array */
  for (i = 0; i < pylen; i++) {
    clist[i] = PyInt_AsLong(PyList_GetItem(pylist, i));
    if (PyErr_Occurred()) {
      return -1;
    }
  }

  return 0;
}

static int prepare_AB(Py_ssize_t *lenA, PyObject *flistA, uint16_t ***fpA, 
		      PyObject *clistA, int **countA, 
		      Py_ssize_t *lenB, PyObject *flistB, uint16_t ***fpB, 
		      PyObject *clistB, int **countB)
{
  Py_ssize_t dblenA, dblenB;
  int dbleni;

  /* Get list lengths */
  *lenA = PyList_Size(flistA);
  *lenB = PyList_Size(flistB);
  if ((*lenA == 0) || (*lenB == 0)) {
    PyErr_SetString(PyExc_ValueError, 
        "Nothing to do since at least one of the lists is empty!");
    return -1;
  }
  if ((*lenA != PyList_Size(clistA)) || (*lenB != PyList_Size(clistB))) {
    PyErr_SetString(PyExc_ValueError, 
        "Bitcount list sizes must match fingerprint list sizes!");
    return -1;
  }

  /* Create c arrays to hold python lists */
  *fpA = (uint16_t **) malloc(*lenA * sizeof(uint16_t *));
  *fpB = (uint16_t **) malloc(*lenB * sizeof(uint16_t *));
  *countA = (int *) malloc(*lenA * sizeof(int));
  *countB = (int *) malloc(*lenB * sizeof(int));
  if ((*fpA == NULL) || (*fpB == NULL) || 
      (*countA == NULL) || (*countB == NULL)) {
    return -2;
  }    

  /* Move listA from python lists to c arrays */  
  dblenA = copy_fp_list(flistA, *lenA, *fpA);
  if (dblenA == -1) {
    return -1;
  }
  if (copy_count_list(clistA, *lenA, *countA) == -1) {
    return -1;
  }

  /* Move listB from python lists to c arrays */  
  dblenB = copy_fp_list(flistB, *lenB, *fpB);
  if (dblenB == -1) {
    return -1;
  }
  else if (dblenA != dblenB) {
    PyErr_SetString(PyExc_ValueError, 
		    "All fingerprints must be the same length!");
    return -1;
  }
  if (copy_count_list(clistB, *lenB, *countB) == -1) {
    return -1;
  }

  dbleni = Py_SAFE_DOWNCAST(dblenA, Py_ssize_t, int);
  return dbleni;
}

/* Calculate number of bits on in raw binary fingerprints
 *
 * Usage: tanimoto.bitcounts(fplist)
 * fplist is a python list of fingerprint raw binary strings
 * returns a python list containing the bitcounts 
 */
static PyObject* tanimoto_bitcounts(PyObject *self, PyObject *args)
{
  PyObject  *flist;
  Py_ssize_t len;
  uint16_t **fp = NULL;
  Py_ssize_t dblen;
  int dbleni;
  Py_ssize_t i;
  int *bcs = NULL;
  PyObject *result;

  /* Parse function argument */
  if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &flist))
    return NULL;

  /* Get list length */
  len = PyList_Size(flist);
  if (len == 0) {
    PyErr_SetString(PyExc_ValueError, 
		    "Nothing to do since the list is empty!");
    return NULL;
  }

  /* Create c array to hold python list */
  fp = (uint16_t **) malloc(len * sizeof(uint16_t *));
  if (fp == NULL) {
    return PyErr_NoMemory();
  }

  /* Move list from python list to c array */
  dblen = copy_fp_list(flist, len, fp);
  if (dblen == -1) {
    free(fp);
    return NULL;
  }
  dbleni = Py_SAFE_DOWNCAST(dblen, Py_ssize_t, int);

  /* Calculate bitcounts into a c array (releasing python GIL) */
  Py_BEGIN_ALLOW_THREADS
  bcs = (int *) malloc(len * sizeof(int));
  if (bcs == NULL) {
    free(fp);
    Py_BLOCK_THREADS
    return PyErr_NoMemory();
  }
  for (i = 0; i < len; i++) {
    bcs[i] = bitcount(fp[i], dbleni);
  }
  Py_END_ALLOW_THREADS

  /* Move bitcounts into python list from c array */
  result = PyList_New(len);
  for (i = 0; i < len; i++) {
    if (PyList_SetItem(result, i, PyInt_FromLong(bcs[i])) == -1) {
      free(fp);
      free(bcs);
      Py_DECREF(result);
      return NULL;
    }
  }

  /* Clean up c arrays and return bitcount list */
  free(fp);
  free(bcs);
  return result;
}

int compare_double(const void *a, const void *b)
{
  double f = *((double *) a);
  double g = *((double *) b);

  if (f > g) return 1;
  if (g > f) return -1;
  return 0;
} 

/* Return the most similar fraction of listB, where similarity is
 *     determined by higher maximum tanimoto to listA
 *
 * Usage: tanimoto.fractionfilterB(fraction, fplistA, countsA, fplistB, countsB)
 * fraction is a double indicating the fraction of listB to return
 * fplistA, fplistB are python lists containing binary fps in python strings
 * countsA, countsB are python lists containing bitcounts for the 
 *     fingerprints in the fplists
 * returns a python list containing indices in listB that have a the highest
 *     maximum tanimotos to listA, to the requested fraction of listB
 */
static PyObject* tanimoto_fractionfilterB(PyObject *self, PyObject *args)
{
  PyObject  *flistA, *clistA, *flistB, *clistB;
  Py_ssize_t lenA, lenB;
  uint16_t **fpA = NULL, **fpB = NULL;
  int *countA = NULL, *countB = NULL;
  int dbleni;
  Py_ssize_t i, j;
  double fraction, threshold;
  double tan, maxtan;
  double *max_tanis, *max_tanis_copy;
  Py_ssize_t *indices;
  Py_ssize_t idx, idlen = 0;
  PyObject *result;

  /* Parse function arguments */
  if (!PyArg_ParseTuple(args, "dO!O!O!O!", &fraction, 
                        &PyList_Type, &flistA, &PyList_Type, &clistA, 
                        &PyList_Type, &flistB, &PyList_Type, &clistB))
    return NULL;

  if ((fraction > 1.0) || (fraction < 0.0)) {
    PyErr_SetString(PyExc_ValueError, 
		    "Fraction not between 0.0 and 1.0.");
    return NULL;
  }

  /* Move fingerprint python lists into c arrays */
  dbleni = prepare_AB(&lenA, flistA, &fpA, clistA, &countA, 
		      &lenB, flistB, &fpB, clistB, &countB);
  if (dbleni == -2) {
    FreeAB
    return PyErr_NoMemory();
  }
  else if (dbleni == -1) {
    FreeAB
    return NULL;
  }

  Py_BEGIN_ALLOW_THREADS
  idx = lenB * (1.0 - fraction);
  max_tanis = (double *) malloc(lenB * sizeof(double));
  if (max_tanis == NULL) {
    FreeAB
    free(max_tanis);
    Py_BLOCK_THREADS
    return PyErr_NoMemory();
  }
  for (j = 0; j < lenB; j++) {
    maxtan = 0.0;
    for (i = 0; i < lenA; i++) {
      tan = tanimoto_coeff(fpA[i], countA[i], fpB[j], countB[j], dbleni);
      if (tan > maxtan) {
        maxtan = tan;
      }
    }
    max_tanis[j] = maxtan;
  }

  max_tanis_copy = (double *) malloc(lenB * sizeof(double));
  if (max_tanis_copy == NULL) {
    FreeAB
    free(max_tanis);
    free(max_tanis_copy);
    Py_BLOCK_THREADS
    return PyErr_NoMemory();
  }
  for (j = 0; j < lenB; j++) {
    max_tanis_copy[j] = max_tanis[j];
  }
  /* Using stdlib qsort for now, even with its function call overhead */
  qsort(max_tanis_copy, lenB, sizeof(double), compare_double);
  threshold = max_tanis_copy[idx];
  free(max_tanis_copy);

  indices = (Py_ssize_t *) malloc(lenB * sizeof(Py_ssize_t));
  if (indices == NULL) {
    FreeAB
    free(max_tanis);
    free(indices);
    Py_BLOCK_THREADS
    return PyErr_NoMemory();
  }
  for (j = 0; j < lenB; j++) {
    if (max_tanis[j] > threshold) {
      indices[idlen++] = j;
    }
  }
  free(max_tanis);
  Py_END_ALLOW_THREADS

  /* Move listB indices from c array to python list */
  result = PyList_New(idlen);
  for (i = 0; i < idlen; i++) {
    if (PyList_SetItem(result, i, PyInt_FromLong(indices[i])) == -1) {
      FreeAB
      free(indices);
      Py_XDECREF(result);
      return NULL;
    }
  }

  /* Clean up c arrays and return list of B indices that are similar */
  FreeAB
  free(indices);
  return result;
}

/* Return indices of listB that are more similar than cutoff to listA
 *
 * Usage: tanimoto.simfilterB(cutoff, fplistA, countsA, fplistB, countsB)
 * fplistA, fplistB are python lists containing binary fps in python strings
 * countsA, countsB are python lists containing bitcounts for the 
 *     fingerprints in the fplists
 * returns a python list containing each index in listB that is more similar 
 *     than the cutoff to any element in listA 
 */
static PyObject* tanimoto_simfilterB(PyObject *self, PyObject *args)
{
  PyObject  *flistA, *clistA, *flistB, *clistB;
  Py_ssize_t lenA, lenB;
  uint16_t **fpA = NULL, **fpB = NULL;
  int *countA = NULL, *countB = NULL;
  int dbleni;
  Py_ssize_t i, j;
  double cutoff;
  double tan;
  Py_ssize_t *indices;
  Py_ssize_t idlen = 0;
  PyObject *result;

  /* Parse function arguments */
  if (!PyArg_ParseTuple(args, "dO!O!O!O!", &cutoff, 
                        &PyList_Type, &flistA, &PyList_Type, &clistA, 
                        &PyList_Type, &flistB, &PyList_Type, &clistB))
    return NULL;

  /* Move fingerprint python lists into c arrays */
  dbleni = prepare_AB(&lenA, flistA, &fpA, clistA, &countA, 
		      &lenB, flistB, &fpB, clistB, &countB);
  if (dbleni == -2) {
    FreeAB
    return PyErr_NoMemory();
  }
  else if (dbleni == -1) {
    FreeAB
    return NULL;
  }

  /* Calculate indices of listB that are more similar than cutoff to listA */
  /*   place results into c array (releasing python GIL) */
  Py_BEGIN_ALLOW_THREADS
  indices = (Py_ssize_t *) malloc(lenB * sizeof(Py_ssize_t));
  if (indices == NULL) {
    FreeAB
    free(indices);
    Py_BLOCK_THREADS
    return PyErr_NoMemory();
  }
  for (j = 0; j < lenB; j++) {
    for (i = 0; i < lenA; i++) {
      tan = tanimoto_coeff(fpA[i], countA[i], fpB[j], countB[j], dbleni);
      if (tan > cutoff) {
	indices[idlen++] = j;
	break;
      }
    }
  }
  Py_END_ALLOW_THREADS

  /* Move listB indices from c array to python list */
  result = PyList_New(idlen);
  for (i = 0; i < idlen; i++) {
    if (PyList_SetItem(result, i, PyInt_FromLong(indices[i])) == -1) {
      FreeAB
      free(indices);
      Py_XDECREF(result);
      return NULL;
    }
  }

  /* Clean up c arrays and return list of B indices that are similar */
  FreeAB
  free(indices);
  return result;
}

/* Return indices of listB that are more similar than each cutoff to listA
 *
 * Usage: tanimoto.batchfilterB(cutlist, fplistA, countsA, fplistB, countsB)
 * cutlist is a python list of integer cutoffs (in the range 0-100)
 * fplistA, fplistB are python lists containing binary fps in python strings
 * fplistA, fplistB are python lists containing binary fps in python strings
 * countsA, countsB are python lists containing bitcounts for the 
 *     fingerprints in the fplists
 * returns a python list of lists where each cutoff in cutlist has a list
 *     of each index in listB that is more similar than that cutoff to any 
 *     element in listA
 */
static PyObject* tanimoto_batchfilterB(PyObject *self, PyObject *args)
{
  PyObject  *cutlist, *flistA, *clistA, *flistB, *clistB;
  Py_ssize_t cutlen, lenA, lenB;
  double *cuts = NULL;
  uint16_t **fpA = NULL, **fpB = NULL;
  int *countA = NULL, *countB = NULL;
  int dbleni;
  Py_ssize_t i, j, k;
  double tan, maxtan;
  Py_ssize_t **indices;
  Py_ssize_t *idlen;
  PyObject *result, *sub;

  /* Parse function arguments */
  if (!PyArg_ParseTuple(args, "O!O!O!O!O!", &PyList_Type, &cutlist, 
                        &PyList_Type, &flistA, &PyList_Type, &clistA, 
                        &PyList_Type, &flistB, &PyList_Type, &clistB))
    return NULL;

  /* Get cutlist length */
  cutlen = PyList_Size(cutlist);
  if (cutlen == 0) {
    PyErr_SetString(PyExc_ValueError, 
        "Nothing to do since the cutoff list is empty!");
    return NULL;
  }
  /* Create c array to hold cutlist */
  cuts = (double *) malloc(cutlen * sizeof(double));
  if (cuts == NULL) {
    return PyErr_NoMemory();
  }    
  /* Move cutlist elements from python list to c array */
  for (i = 0; i < cutlen; i++) {
    cuts[i] = PyFloat_AsDouble(PyList_GetItem(cutlist, i));
    if (PyErr_Occurred()) {
      free(cuts);
      return NULL;
    }
  }

  /* Move fingerprint python lists into c arrays */
  dbleni = prepare_AB(&lenA, flistA, &fpA, clistA, &countA, 
		      &lenB, flistB, &fpB, clistB, &countB);
  if (dbleni == -2) {
    free(cuts);
    FreeAB
    return PyErr_NoMemory();
  }
  else if (dbleni == -1) {
    free(cuts);
    FreeAB
    return NULL;
  }

  /* Generate c 2-dimensional array to hold filtered indices for each cutoff */
  Py_BEGIN_ALLOW_THREADS
  idlen = (Py_ssize_t *) malloc(cutlen * sizeof(Py_ssize_t));
  indices = (Py_ssize_t **) malloc(cutlen * sizeof(Py_ssize_t *));
  if ((idlen == NULL) || (indices == NULL)) {
    free(cuts);
    FreeAB
    free(idlen);
    free(indices);
    Py_BLOCK_THREADS
    return PyErr_NoMemory();
  }
  for (i = 0; i < cutlen; i++) {
    idlen[i] = 0;
    indices[i] = (Py_ssize_t *) malloc(lenB * sizeof(Py_ssize_t));
    if (indices[i] == NULL) {
      free(cuts);
      FreeAB
      free(idlen);
      for (j = 0; j < i; j++) {
	free(indices[j]);
      }
      free(indices);
      Py_BLOCK_THREADS
      return PyErr_NoMemory();
    }
  }

  /* Calculate filtered indices into c 2-dimensional array */
  for (j = 0; j < lenB; j++) {
    maxtan = 0.0;
    for (i = 0; i < lenA; i++) {
      tan = tanimoto_coeff(fpA[i], countA[i], fpB[j], countB[j], dbleni);
      if (tan > maxtan) {
	  maxtan = tan;
      }
    }
    for (k = 0; k < cutlen; k++) {
      if (maxtan > cuts[k]) {
	indices[k][idlen[k]] = j;
	idlen[k] += 1;
      }
    }
  }

  Py_END_ALLOW_THREADS

  /* Move calculated indices from c 2-d array to python list of lists */
  result = PyList_New(cutlen);
  for (i = 0; i < cutlen; i++) {
    sub = PyList_New(idlen[i]);
    for (j = 0; j < idlen[i]; j++) {
      if (PyList_SetItem(sub, j, PyInt_FromLong(indices[i][j])) == -1) {
	free(cuts);
	FreeAB
	free(idlen);
	for (k = 0; k < cutlen; k++) {
	  free(indices[k]);
	}
	free(indices);
	Py_DECREF(sub);
	Py_DECREF(result);
	return NULL;
      }
    }
    if (PyList_SetItem(result, i, sub) == -1) {
      free(cuts);
      FreeAB
      free(idlen);
      for (k = 0; k < cutlen; k++) {
	free(indices[k]);
      }
      free(indices);
      Py_XDECREF(sub);
      Py_XDECREF(result);
      return NULL;
    }
  }

  /* Clean up c arrays and return list of B indices that are similar */
  free(cuts);
  FreeAB
  free(idlen);
  for (k = 0; k < cutlen; k++) {
    free(indices[k]);
  }
  free(indices);
  return result;
}

/* Calculate the pairwise tanimoto matrix for fplistA versus fplistB
 * 
 * Usage: tanimoto.matrix(fplistA, countsA, fplistB, countsB)
 * fplistA, fplistB are python lists containing binary fps in python strings
 * countsA, countsB are python lists containing bitcounts for the 
 *     fingerprints in the fplists
 * returns a python list of lists containing the tanimoto matrix of all
 *     pairwise comparisions of fplistA elements with fplistB elements
 */
static PyObject* tanimoto_matrix(PyObject *self, PyObject *args)
{
  PyObject  *flistA, *clistA, *flistB, *clistB;
  Py_ssize_t lenA, lenB;
  uint16_t **fpA = NULL, **fpB = NULL;
  int *countA = NULL, *countB = NULL;
  int dbleni;
  Py_ssize_t i, j, k;
  double **tm;
  PyObject *result, *sub;


  /* Parse function arguments */
  if (!PyArg_ParseTuple(args, "O!O!O!O!", 
                        &PyList_Type, &flistA, &PyList_Type, &clistA, 
                        &PyList_Type, &flistB, &PyList_Type, &clistB))
    return NULL;

  /* Move fingerprint python lists into c arrays */
  dbleni = prepare_AB(&lenA, flistA, &fpA, clistA, &countA, 
		      &lenB, flistB, &fpB, clistB, &countB);
  if (dbleni == -2) {
    FreeAB
    return PyErr_NoMemory();
  }
  else if (dbleni == -1) {
    FreeAB
    return NULL;
  }

  /* Calculate pairwise tanimoto matrix into c 2-dimensional array */
  Py_BEGIN_ALLOW_THREADS
  tm = (double **) malloc(lenA * sizeof(double *));
  if (tm == NULL) {
    FreeAB
    free(tm);
    Py_BLOCK_THREADS
    return PyErr_NoMemory();
  }
  for (i = 0; i < lenA; i++) {
    tm[i] = (double *) malloc(lenB * sizeof(double));
    if (tm[i] == NULL) {
      FreeAB
      for (j = 0; j < i; j++) {
	free(tm[j]);
      }
      free(tm);
      Py_BLOCK_THREADS
      return PyErr_NoMemory();
    }
    for (j = 0; j < lenB; j++) {
      tm[i][j] = tanimoto_coeff(fpA[i], countA[i], 
				fpB[j], countB[j], dbleni);
    }
  }
  Py_END_ALLOW_THREADS

  /* Move c 2-dimensional array into python list of lists */
  result = PyList_New(lenA);
  for (i = 0; i < lenA; i++) {
    sub = PyList_New(lenB);
    for (j = 0; j < lenB; j++) {
      if (PyList_SetItem(sub, j, PyFloat_FromDouble(tm[i][j])) == -1) {
	FreeAB
	for (k = 0; k < lenA; k++) {
	  free(tm[k]);
	}
	free(tm);
        Py_XDECREF(sub);
        Py_XDECREF(result);
        return NULL;
      }
    }
    if (PyList_SetItem(result, i, sub) == -1) {
      FreeAB
      for (k = 0; k < lenA; k++) {
	free(tm[k]);
      }
      free(tm);
      Py_XDECREF(sub);
      Py_XDECREF(result);
      return NULL;
    }
  }

  /* Clean up c arrays and return tanimoto maxtrix */
  FreeAB
  for (k = 0; k < lenA; k++) {
    free(tm[k]);
  }
  free(tm);
  return result;
}

/* Load ascii fp list and return raw fp list and bitcounts
 *
 * Usage: tanimoto.loadascii(alist)
 * alist is a python list of the input Daylight ASCII fingerprints
 * returns a tuple of two lists: (rlist, bclist)
 * rlist is a python list containing the converted raw binary strings
 * bclist is a python list containing the bitcounts
 */
static PyObject* tanimoto_loadascii(PyObject *self, PyObject *args)
{
  PyObject  *alist;
  Py_ssize_t len;
  char **asc;
  Py_ssize_t i, j;
  int dblen;
  char **rstr;
  int *rlen, *bcs;
  PyObject *rlist, *bclist, *x, *ret;

  /* Parse function argument */
  if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &alist))
    return NULL;

  /* Get list length */
  len = PyList_Size(alist);
  if (len == 0) {
    PyErr_SetString(PyExc_ValueError, 
		    "Nothing to do since the list is empty!");
    return NULL;
  }

  /* Create c array to hold python list */
  asc = (char **) malloc(len * sizeof(char *));
  if (asc == NULL) {
    return PyErr_NoMemory();
  }    

  /* Move list elements from python list to c array */
  for (i = 0; i < len; i++) {
    asc[i] = PyString_AsString(PyList_GetItem(alist, i));
    if ( asc[i] == NULL) {
      free(asc);
      return NULL;
    }
  }

  /* Convert the strings and calculate bitcount in c arrays */
  Py_BEGIN_ALLOW_THREADS
  rstr = (char **) malloc(len * sizeof(char *));
  rlen = (int *) malloc(len * sizeof(int));
  bcs = (int *) malloc(len * sizeof(int));
  if ((rstr == NULL) || (rlen == NULL) || (bcs == NULL)) {
    free(asc);
    free(rstr);
    free(rlen);
    free(bcs);
    Py_BLOCK_THREADS
    return PyErr_NoMemory();
  }

  rstr[0] = ascii2binary(asc[0], &(rlen[0]));
  if (rstr[0] == NULL) {
      free(asc);
      free(rstr);
      free(rlen);
      free(bcs);
      Py_BLOCK_THREADS
      PyErr_SetString(PyExc_ValueError, 
		      "Cannot convert first ascii string!");
      return NULL;
  }
  dblen = rlen[0]/2;

  for (i = 1; i < len; i++) {
    rstr[i] = ascii2binary(asc[i], &(rlen[i]));
    if (rstr[i] == NULL) {
      free(asc);
      for (j = 0; j < i; j++) {
	free(rstr[j]);
      }
      free(rstr);
      free(rlen);
      free(bcs);
      Py_BLOCK_THREADS
      PyErr_SetString(PyExc_ValueError, 
		      "Cannot convert ascii strings to raw binary strings!");
      return NULL;
    }
    if ((dblen*2) != rlen[i]) {
      free(asc);
      for (j = 0; j < i; j++) {
	free(rstr[j]);
      }
      free(rstr);
      free(rlen);
      free(bcs);
      Py_BLOCK_THREADS
      PyErr_SetString(PyExc_ValueError, 
		      "All fingerprints must be the same length!");
      return NULL;
    }
  }

  for (i = 0; i < len; i++) {
    bcs[i] = bitcount((uint16_t *) rstr[i], dblen);
  }
  Py_END_ALLOW_THREADS

  /* Move raw strings into python list from c arrays */
  rlist = PyList_New(len);
  for (i = 0; i < len; i++) {
    x = PyString_FromStringAndSize(rstr[i], (Py_ssize_t) rlen[i]);
    if (PyList_SetItem(rlist, i, x) == -1) {
      free(asc);
      for (j = 0; j < len; j++) {
	free(rstr[j]);
      }
      free(rstr);
      free(rlen);
      free(bcs);
      Py_XDECREF(rlist);
      return NULL;
    }
  }

  /* Move bitcounts into python list from c array */
  bclist = PyList_New(len);
  for (i = 0; i < len; i++) {
    if (PyList_SetItem(bclist, i, PyInt_FromLong(bcs[i])) == -1) {
      free(asc);
      for (j = 0; j < len; j++) {
	free(rstr[j]);
      }
      free(rstr);
      free(rlen);
      free(bcs);
      Py_XDECREF(rlist);
      Py_XDECREF(bclist);
      return NULL;
    }
  }

  /* Clean up and return the lists of raw fingerprints and bitcounts */
  free(asc);
  for (i = 0; i < len; i++) {
    free(rstr[i]);
  }
  free(rstr);
  free(rlen);
  free(bcs);
  ret = Py_BuildValue("(OO)", rlist, bclist);
  Py_XDECREF(rlist);
  Py_XDECREF(bclist);
  return ret;
}

/* Convert Daylight ASCII fingerprints to raw binary fingerprints
 *
 * Usage: tanimoto.asc2bin(alist)
 * alist is a python list of the input Daylight ASCII fingerprints
 * returns a python list containing the converted raw binary strings
 */
static PyObject* tanimoto_asc2bin(PyObject *self, PyObject *args)
{
  PyObject  *alist;
  Py_ssize_t len;
  char **asc;
  Py_ssize_t i, j;
  char **rstr;
  int *rlen;
  PyObject *result, *x;

  /* Parse function argument */
  if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &alist))
    return NULL;

  /* Get list length */
  len = PyList_Size(alist);
  if (len == 0) {
    PyErr_SetString(PyExc_ValueError, 
		    "Nothing to do since the list is empty!");
    return NULL;
  }

  /* Create c array to hold python list */
  asc = (char **) malloc(len * sizeof(char *));
  if (asc == NULL) {
    return PyErr_NoMemory();
  }    

  /* Move list elements from python list to c array */
  for (i = 0; i < len; i++) {
    asc[i] = PyString_AsString(PyList_GetItem(alist, i));
    if (asc[i] == NULL) {
      free(asc);
      return NULL;
    }
  }

  /* Convert the strings in c arrays (releasing python GIL) */
  Py_BEGIN_ALLOW_THREADS
  rstr = (char **) malloc(len * sizeof(char *));
  rlen = (int *) malloc(len * sizeof(int));
  if ((rstr == NULL) || (rlen == NULL)) {
    free(asc);
    free(rstr);
    free(rlen);
    Py_BLOCK_THREADS
    return PyErr_NoMemory();
  }
  for (i = 0; i < len; i++) {
    rstr[i] = ascii2binary(asc[i], &(rlen[i]));
  }
  Py_END_ALLOW_THREADS

  /* Move raw strings into python list from c arrays */
  result = PyList_New(len);
  for (i = 0; i < len; i++) {
    x = PyString_FromStringAndSize(rstr[i], (Py_ssize_t) rlen[i]);
    if (PyList_SetItem(result, i, x) < 0) {
      free(asc);
      for (j = 0; j < len; j++) {
	free(rstr[j]);
      }
      free(rstr);
      free(rlen);
      Py_XDECREF(result);
      return NULL;
    }
  }

  /* Clean up c arrays and return the list of raw fingerprints */
  free(asc);
  for (i = 0; i < len; i++) {
    free(rstr[i]);
  }
  free(rstr);
  free(rlen);
  return result;
}

/* Convert raw binary fingerprints to Daylight ASCII fingerprints 
 *
 * Usage: tanimoto.bin2asc(blist)
 * blist is a python list of the input raw binary strings
 * returns a python list containing the converted Daylight ASCII fingerprints
 */
static PyObject* tanimoto_bin2asc(PyObject *self, PyObject *args)
{
  PyObject  *blist;
  Py_ssize_t len;
  char **rstr;
  Py_ssize_t *rlen;
  Py_ssize_t i, j;
  int rleni;
  char **astr;
  PyObject *result, *x;

  /* Parse function argument */
  if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &blist))
    return NULL;

  /* Get list length */
  len = PyList_Size(blist);
  if (len == 0) {
    PyErr_SetString(PyExc_ValueError, 
		    "Nothing to do since the list is empty!");
    return NULL;
  }

  /* Create c array to hold python list */
  rstr = (char **) malloc(len * sizeof(char *));
  rlen = (Py_ssize_t *) malloc(len * sizeof(Py_ssize_t));
  if ((rstr == NULL) || (rlen == NULL)) {
    free(rstr);
    free(rlen);
    return PyErr_NoMemory();
  }    

  /* Move list elements from python list to c array */
  for (i = 0; i < len; i++) {
    x = PyList_GetItem(blist, i);
    if (PyString_AsStringAndSize(x, &(rstr[i]), &(rlen[i])) == -1 ) {
      free(rstr);
      free(rlen);
      return NULL;
    }
  }

  /* Convert the strings in c arrays (releasing python GIL) */
  Py_BEGIN_ALLOW_THREADS
  astr = (char **) malloc(len * sizeof(char *));
  if (astr == NULL) {
    free(rstr);
    free(rlen);
    Py_BLOCK_THREADS
    return PyErr_NoMemory();
  }
  for (i = 0; i < len; i++) {
    rleni = Py_SAFE_DOWNCAST(rlen[i], Py_ssize_t, int);
    astr[i] = binary2ascii(rstr[i], rleni);
  }
  Py_END_ALLOW_THREADS

  /* Move raw strings into python list from c arrays */
  result = PyList_New(len);
  for (i = 0; i < len; i++) {
    if (PyList_SetItem(result, i, PyString_FromString(astr[i])) == -1) {
      free(rstr);
      free(rlen);
      for (j = 0; j < len; j++) {
	free(astr[j]);
      }
      free(astr);
      Py_XDECREF(result);
      return NULL;
    }
  }

  /* Clean up c arrays and return the list of raw fingerprints */
  free(rstr);
  free(rlen);
  for (i = 0; i < len; i++) {
    free(astr[i]);
  }
  free(astr);
  return result;
}

/* PyMethodDef structure relating python functions with c functions */
static PyMethodDef tanimoto_methods[] = {
    {"matrix", tanimoto_matrix, METH_VARARGS,
     "Calculate matrix of tanimoto coefficients."}, 
    {"fractionfilterB", tanimoto_fractionfilterB, METH_VARARGS,
     "Return top fraction of listB indices ordered by decreasing maximum tanimoto to listA"}, 
    {"simfilterB", tanimoto_simfilterB, METH_VARARGS,
     "Return indices of listB that are more similar than cutoff to listA"}, 
    {"batchfilterB", tanimoto_batchfilterB, METH_VARARGS,
     "Return indices of listB that are more similar than each cutoff to listA"}, 
    {"loadascii", tanimoto_loadascii, METH_VARARGS,
     "Load ascii fp list and return raw fp list and bitcounts"}, 
    {"bitcounts", tanimoto_bitcounts, METH_VARARGS,
     "Calculate number of bits on in raw binary fingerprints"}, 
    {"asc2bin", tanimoto_asc2bin, METH_VARARGS,
     "Convert Daylight ASCII fingerprints to raw binary fingerprints"}, 
    {"bin2asc", tanimoto_bin2asc, METH_VARARGS, 
     "Convert raw binary fingerprints to Daylight ASCII fingerprints"}, 
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

/* Initialize tanimoto module with the PyMethodDef structure 
 * Then call init_BitCountArray() to initialize fast tanimoto lookup table
 */
PyMODINIT_FUNC inittanimoto(void)
{
  /* tanimotoX numbering is a debugging hack, remove later */
  Py_InitModule("tanimoto", tanimoto_methods);
  init_BitCountArray();
}
