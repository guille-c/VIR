#ifndef VIS_LISTA
#define VIS_LISTA

#include "estructuras_malla.h"
#include "nodoLista.h"

#ifdef __cplusplus
extern "C" {
#endif  

struct lista {
  NodoLista *primero;
};

Lista *newLista();
void insertarNodo(Lista *l, NodoLista *n);
void eliminarNodo(Lista *l, NodoLista *n);
void eliminarLista(Lista *l);
void imprimirListaInt(Lista *l);

#ifdef __cplusplus
}
#endif  


#endif
