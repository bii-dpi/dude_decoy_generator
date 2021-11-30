/* fastdl: Module for FAST daylight toolkit use from python
 * 
 * Function List:
 *
 * fingerprints([(key1, smi1), ...])
 * Calculate ascii fingerprints from smiles 
 *
 * Revisions:
 * Mysinger 8/06 Created
 */

#include <Python.h>
#include <dt_smiles.h>
#include <dt_finger.h>

#if PY_VERSION_HEX < 0x02050000 && !defined(PY_SSIZE_T_MIN)
typedef int Py_ssize_t;
#define PY_SSIZE_T_MAX INT_MAX
#define PY_SSIZE_T_MIN INT_MIN
#endif

/* Calculate ascii fingerprints from smiles 
 *
 * Usage: fastdl.fingerprints(input_list)
 * input_list = [(key_object1, smiles1), ...] 
 * return_list = [(key_object1, ascii_fp1), ...]
 *
 * Requires the Daylight Toolkit be installed and DY_ROOT is set
 */
static PyObject* fastdl_fingerprints(PyObject *self, PyObject *args)
{
  PyObject *smilist, *t;
  Py_ssize_t i, j, slen;
  PyObject **ids, **refs;
  char **smiles;
  dt_Handle mol, fp, dth;
  char *rstr, *tempstr;
  int rlen, templen;
  char **astr;
  PyObject *result;

  /* Parse function argument */
  if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &smilist))
    return NULL;
  slen = PyList_Size(smilist);
  if (slen == 0) {
    PyErr_SetString(PyExc_ValueError, 
      "Nothing to do since the list is empty!");
    return NULL;
  }

  /* Create c arrays to hold python list */
  refs = (PyObject **) malloc(slen * sizeof(PyObject *));
  ids = (PyObject **) malloc(slen * sizeof(PyObject *));
  smiles = (char **) malloc(slen * sizeof(char *));
  if ((refs == NULL) || (ids == NULL) || (smiles == NULL)) {
    free(ids);
    free(smiles);
    free(refs);
    return PyErr_NoMemory();
  }

  /* Move elements from python list of tuples to c arrays */
  for (i = 0; i < slen; i++) {
    refs[i] = PySequence_Tuple(PyList_GetItem(smilist, i));
    /*  if (!PyArg_ParseTuple(refs[i], "Os", &(ids[i]), &(smiles[i]))) { */
    if (!PyArg_ParseTuple(refs[i], "sO", &(smiles[i]), &(ids[i]))) {
      free(ids);
      free(smiles);
      for (j = 0; j < i; j++) {
	Py_DECREF(refs[j]);
      }
      Py_XDECREF(refs[i]);
      free(refs);
      return NULL;
    }
  }

  /* Use Daylight Toolkit to calculate raw and ascii fingerprint strings */
  Py_BEGIN_ALLOW_THREADS
  astr = (char **) malloc(slen * sizeof(char *));
  if (astr == NULL) {
    free(ids);
    free(smiles);
    Py_BLOCK_THREADS
    for (j = 0; j < slen; j++) {
      Py_DECREF(refs[j]);
    }
    free(refs);
    return PyErr_NoMemory();
  }
  for (i = 0; i < slen; i++) {
    mol = dt_smilin(strlen(smiles[i]), smiles[i]);
    if (mol == NULL_OB) {
      printf("Error generating molecule on smiles \"%s\", skipping!\n", 
	     smiles[i]);
      astr[i] = NULL;
      continue;
    }
    fp = dt_fp_generatefp(mol, 0, 7, 2048);
    if (fp == NULL_OB) {
      printf("Error generating fingerprint on smiles \"%s\", skipping!\n", 
	     smiles[i]);
      astr[i] = NULL;
      continue;
    }
    rstr = dt_stringvalue(&rlen, fp);
    if (rstr == NULL) {
      printf("Error getting raw string on smiles \"%s\", skipping!\n", 
	     smiles[i]);
      astr[i] = NULL;
      continue;
    }
    dth = dt_binary2ascii(rlen, rstr);
    if (dth == NULL_OB) {
      printf("Error generating ascii string on smiles \"%s\", skipping!\n", 
	     smiles[i]);
      astr[i] = NULL;
      continue;
    }
    tempstr = dt_stringvalue(&templen, dth);
    if (tempstr == NULL) {
      printf("Error getting ascii string on smiles \"%s\", skipping!\n", 
	     smiles[i]);
      astr[i] = NULL;
      continue;
    }
    astr[i] = (char *) malloc((templen + 1) * sizeof(char));
    if (astr[i] == NULL) {
      dt_dealloc(mol);
      dt_dealloc(fp);
      dt_dealloc(dth);
      free(ids);
      free(smiles);
      for (j = 0; j < i; j++) {
	free(astr[j]);
      }
      free(astr);
      Py_BLOCK_THREADS
      for (j = 0; j < slen; j++) {
        Py_DECREF(refs[j]);
      }
      free(refs);
      return PyErr_NoMemory();
    }
    memcpy(astr[i], tempstr, templen);
    astr[i][templen] = '\0';
    dt_dealloc(mol);
    dt_dealloc(fp);
    dt_dealloc(dth);
  }
  Py_END_ALLOW_THREADS

  /* Move results to a list of python tuples from c arrays */
  result = PyList_New(0);
  for (i = 0; i < slen ; i++) {
    /* t = Py_BuildValue("(Os)", ids[i], astr[i]); */
    t = Py_BuildValue("(sO)", astr[i], ids[i]);
    if (PyList_Append(result, t) < 0) {
      Py_XDECREF(t);
      Py_DECREF(result);
      free(ids);
      free(smiles);
      for (i = 0; i < slen; i++) {
	free(astr[i]);
      }
      free(astr);
      for (j = 0; j < slen; j++) {
        Py_DECREF(refs[j]);
      }
      free(refs);
      return NULL;
    }
    Py_DECREF(t);
  }

  /* Clean up c arrays and return fingerprint list */
  free(ids);
  free(smiles);
  for (i = 0; i < slen; i++) {
    free(astr[i]);
  }
  free(astr);
  for (j = 0; j < slen; j++) {
    Py_DECREF(refs[j]);
  }
  free(refs);
  return result;
}

/* PyMethodDef structure relating python functions with c functions */
static PyMethodDef fastdl_methods[] = {
    {"fingerprints", fastdl_fingerprints, METH_VARARGS,
     "Calculate fingerprint strings given ids and smiles."}, 
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

/* Initialize fastdl module with the PyMethodDef structure */
PyMODINIT_FUNC initfastdl(void)
{
  Py_InitModule("fastdl", fastdl_methods);
}
