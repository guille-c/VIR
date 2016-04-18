#include "lista.h"
#include <stdio.h>
#include "punto.h"
#include <malloc.h>

Lista *newLista() {
  Lista *l =
    (Lista *) malloc(sizeof(Lista));
  l->primero = NULL;
  return l;
}

void insertarNodo(Lista *l, NodoLista *n) {
  if(l->primero == NULL) {
    l->primero = n;
    return;
  }
  else {
    n->sgte = l->primero;
    l->primero->ante = n;
    l->primero = n;
  }
}

void eliminarNodo(Lista *l, NodoLista *n) {
  if (n == NULL) {
    return;
  }
  if (n == l->primero) {
    l->primero = n->sgte;
  }
  if(n->sgte != NULL) {
    n->sgte->ante = n->ante;
  }
  if(n->ante != NULL) {
    n->ante->sgte = n->sgte;
  }
  eliminarNodoLista(n);
}

void eliminarLista(Lista *l) {
  while(l->primero != NULL) {
    eliminarNodo(l, l->primero);
  }
  free(l);
}

void imprimirListaInt(Lista *l) {
  NodoLista *n = l->primero;

  printf("{ ");
  while (n != NULL) {
    printf("%d ", (int) (n->info));
    n = n->sgte;
  }
  printf("}\n");
}
