#include "nodoLista.h"
#include "punto.h"
#include <malloc.h>

NodoLista *newNodoLista(void *info,
			NodoLista *ante,
			NodoLista *sgte) {
  NodoLista *n = (NodoLista *) malloc(sizeof(NodoLista));
  n->info = info;
  n->ante = ante;
  n->sgte = sgte;
  return n;
}

void eliminarNodoLista(NodoLista *n) {
  //free(n->info);
  if (n != NULL) {
    free(n);
  }
}
