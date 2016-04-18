#ifndef VIS_NODO_LISTA
#define VIS_NODO_LISTA

#include "estructuras_malla.h"

#ifdef __cplusplus
extern "C" {
#endif  

struct nodoLista {
  void *info;
  NodoLista *ante, *sgte;
};

NodoLista *newNodoLista(void *info,
			NodoLista *ante,
			NodoLista *sgte);
void eliminarNodoLista(NodoLista *n);

#ifdef __cplusplus
}
#endif

#endif
