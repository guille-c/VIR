#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mallaVoronoi.h"
#include "rasterizador.h"

#define DEBUG_MALLA 0

MallaVoronoi *newMallaVoronoi() {
  PoligonoVoronoi *pol1, *pol2, *pol3;
  PuntoVoronoi *p, *p1, *p2, *p3;
  AristaVoronoi *a1, *a2, *a3, *aIm1, *aIm2, *aIm3;
  double x, y;
  MallaVoronoi *m = 
    (MallaVoronoi *) malloc(sizeof(MallaVoronoi));
  NodoLista *n;
  
  m->aristas = newLista();
  m->poligonos = newLista();
  m->puntos = newLista();
  m->nPols = 3;

  //pol1 = newPoligonoVoronoi(0.5, 3.0 * sqrt(2.0)/4.0 + 1.0, NULL_PIX, 0);
  pol1 = newPoligonoVoronoi(0.5, 3.0 * sqrt(2.0)/4.0 + 1.0, 0, 0);
  n = newNodoLista (pol1, NULL, NULL);
  pol1->n = n;
  insertarNodo (m->poligonos, n);

  //pol2 = newPoligonoVoronoi(-3.0 * sqrt(6.0)/4.0 + 0.5, 
  //-3.0 * sqrt(2.0)/4.0 + 0.5, NULL_PIX, 1);
  pol2 = newPoligonoVoronoi(-3.0 * sqrt(6.0)/4.0 + 0.5, 
			    -3.0 * sqrt(2.0)/4.0 + 0.5, 0, 1);
  n = newNodoLista (pol2, NULL, NULL);
  pol2->n = n;
  insertarNodo (m->poligonos, n);

  //pol3 = newPoligonoVoronoi(3.0 * sqrt(6.0)/4.0 + 0.5, 
  //-3.0 * sqrt(2.0)/4.0 + 0.5, NULL_PIX, 2);
  pol3 = newPoligonoVoronoi(3.0 * sqrt(6.0)/4.0 + 0.5, 
			    -3.0 * sqrt(2.0)/4.0 + 0.5, 0, 2);
  n = newNodoLista (pol3, NULL, NULL);
  pol3->n = n;
  insertarNodo (m->poligonos, n);
  
  /*  pol1 = newPoligonoVoronoi(1.5, 2.6, 1.0);
      pol2 = newPoligonoVoronoi(0.0, 0.5, 2.0);
      pol3 = newPoligonoVoronoi(1.7, 0.0, 3.0);
  */
  p = circuncentro(pol1, pol2, pol3);
  n = newNodoLista (p, NULL, NULL);
  p->n = n;
  insertarNodo (m->puntos, n);
  //printf("p: (%g, %g)\n", p->x - 0.5, p->y);
  
  y = -10.0;
  x = - (pol2->y - pol3->y) / (pol2->x - pol3->x) * (y - p->y) + p->x;
  p1 = newPuntoVoronoi(x, y);
  n = newNodoLista (p1, NULL, NULL);
  p1->n = n;
  insertarNodo (m->puntos, n);

  x = 10.0;
  y = - (pol1->x - pol3->x) / (pol1->y - pol3->y) * (x - p->x) + p->y;
  p2 = newPuntoVoronoi(x, y);
  n = newNodoLista (p2, NULL, NULL);
  p2->n = n;
  insertarNodo (m->puntos, n);

  x = -10.0;
  y = - (pol1->x - pol2->x) / (pol1->y - pol2->y) * (x - p->x) + p->y;
  p3 = newPuntoVoronoi(x, y);
  n = newNodoLista (p3, NULL, NULL);
  p3->n = n;
  insertarNodo (m->puntos, n);

  a1 = newAristaVoronoi (p, p1);
  a2 = newAristaVoronoi (p, p2);
  a3 = newAristaVoronoi (p, p3);
  aIm1 = newAristaVoronoi (p2, p3);
  aIm2 = newAristaVoronoi (p3, p1);
  aIm3 = newAristaVoronoi (p1, p2);
  
  a1->poliDer = pol2;
  a1->poliIzq = pol3;
  a1->cwPred = a3;
  a1->ccwPred = a2;
  a1->cwSucc = aIm3;
  a1->ccwSucc = aIm2;
  calcularSinCos(a1);
  
  a2->poliDer = pol3;
  a2->poliIzq = pol1;
  a2->cwPred = a1;
  a2->ccwPred = a3;
  a2->cwSucc = aIm1;
  a2->ccwSucc = aIm3;
  calcularSinCos(a2);
  
  a3->poliDer = pol1;
  a3->poliIzq = pol2;
  a3->cwPred = a2;
  a3->ccwPred = a1;
  a3->cwSucc = aIm2;
  a3->ccwSucc = aIm1;
  calcularSinCos(a3);

  aIm1->poliDer = NULL;
  aIm1->poliIzq = pol1;
  aIm1->cwPred = aIm3;
  aIm1->ccwPred = a2;
  aIm1->cwSucc = a3;
  aIm1->ccwSucc = aIm2;
  calcularSinCos(aIm1);

  aIm2->poliDer = NULL;
  aIm2->poliIzq = pol2;
  aIm2->cwPred = aIm1;
  aIm2->ccwPred = a3;
  aIm2->cwSucc = a1;
  aIm2->ccwSucc = aIm3;
  calcularSinCos(aIm2);

  aIm3->poliDer = NULL;
  aIm3->poliIzq = pol3;
  aIm3->cwPred = aIm2;
  aIm3->ccwPred = a1;
  aIm3->cwSucc = a2;
  aIm3->ccwSucc = aIm1;
  calcularSinCos(aIm3);

  
  //imprimirArista(a1);
  //imprimirArista(a2);
  //imprimirArista(a3);
  
  p->a = a1;
  p1->a = a1;
  p2->a = a2;
  p3->a = a3;
  
  pol1->a = a2;
  pol2->a = a3;
  pol3->a = a1;
  
  n = newNodoLista(a1, NULL, NULL);
  a1->n = n;
  insertarNodo(m->aristas, n);

  n = newNodoLista(a2, NULL, NULL);
  a2->n = n;
  insertarNodo(m->aristas, n);

  n = newNodoLista(a3, NULL, NULL);
  a3->n = n;
  insertarNodo(m->aristas, n);

  n = newNodoLista(aIm1, NULL, NULL);
  aIm1->n = n;
  insertarNodo(m->aristas, n);

  n = newNodoLista(aIm2, NULL, NULL);
  aIm2->n = n;
  insertarNodo(m->aristas, n);

  n = newNodoLista(aIm3, NULL, NULL);
  aIm3->n = n;
  insertarNodo(m->aristas, n);

  return m;
}

PuntoVoronoi *circuncentro(PoligonoVoronoi *a,
				  PoligonoVoronoi *b,
				  PoligonoVoronoi *c) {
  double xba, yba, xca, yca;
  double balength, calength;
  double zcrossbc;
  double denominator;
  double xcirca, ycirca;
  
  /* Use coordinates relative to point `a' of the triangle. */
  xba = b->x - a->x;
  yba = b->y - a->y;
  xca = c->x - a->x;
  yca = c->y - a->y;
  
  /* Squares of lengths of the edges incident to `a'. */
  balength = xba * xba + yba * yba;
  calength = xca * xca + yca * yca;
  
  /* Cross product of these edges. */
  zcrossbc = xba * yca - xca * yba;
  
  /* Calculate the denominator of the formulae. */
  denominator = 0.5 / (zcrossbc * zcrossbc);
  
  /* Calculate offset (from `a') of circumcenter. */
  xcirca = ((balength * yca - calength * yba) * zcrossbc) * denominator;
  ycirca = (-(balength * xca - calength * xba) * zcrossbc) * denominator;

  return newPuntoVoronoi(a->x + xcirca, a->y + ycirca);
}

void eliminarMallaVoronoi(MallaVoronoi *m) {
  NodoLista *n;
  AristaVoronoi *a;
  PuntoVoronoi *p;
  PoligonoVoronoi *pol;

  if (m->aristas == NULL) {
    return;
  }
  for (n = m->aristas->primero; n != NULL; n = n->sgte) {
    a = (AristaVoronoi *) n->info;
    eliminarNodo (m->aristas, a->n);
    if (a != NULL) {
      eliminarAristaVoronoi (a);
    }
  }
  eliminarLista(m->aristas);
  m->aristas = NULL;

  for (n = m->puntos->primero; n != NULL; n = n->sgte) {
    p = (PuntoVoronoi *) n->info;
    eliminarNodo (m->puntos, p->n);
    if (p != NULL) {
      free (p);
    }
  }
  eliminarLista(m->puntos);

  for (n = m->poligonos->primero; n != NULL; n = n->sgte) {
    pol = (PoligonoVoronoi *) n->info;
    eliminarNodo (m->poligonos, pol->n);
    if (pol != NULL) {
      free (pol);
    }
  }
  eliminarLista(m->poligonos);
  free(m);
}

// Metodo que retorna el poligono mas cercano a (x, y)
PoligonoVoronoi *encontrarPoligono(MallaVoronoi *m, 
				   double x, double y) {
  AristaVoronoi *a;
  
  if (m->aristas == NULL || m->aristas->primero == NULL) {
    return NULL;
  }
  a = m->aristas->primero->info;
  if (a->poliDer != NULL) {
    return encontrarPoligono2(a->poliDer, x, y);
  }
  else if (a->poliIzq != NULL){
    return encontrarPoligono2(a->poliIzq, x, y);
  } 
  return NULL; 
}

/* Metodo que retorna el poligono mas cercano a (x, y) 
 * comenzando en polIni 
 */ 
PoligonoVoronoi *encontrarPoligono2(PoligonoVoronoi *polIni, 
					   double x, double y) { 
  
  PoligonoVoronoi *escogido; 
  AristaVoronoi *a; 
  double d1, d2; 
  
  escogido = NULL;
  //d2 = pow(x - polIni->x, 2) + pow(y - polIni->y, 2);
  d2 = (x - polIni->x) * (x - polIni->x) + (y - polIni->y) * (y - polIni->y);
  a = polIni->a;
  
  //printf ("Buscando en poligono: (%g, %g), x = %g, y = %g, d2 = %g\n", 
  //polIni->x, polIni->y, x, y, d2);
  do {
    //imprimirArista(a);
    if (a->poliDer == polIni) {
      if (a->poliIzq != NULL && 
	  (d1 = ((x - a->poliIzq->x) * (x - a->poliIzq->x) + 
		 (y - a->poliIzq->y) * (y - a->poliIzq->y))) < d2) {
	d2 = d1;
	escogido = a->poliIzq;
      }
      a = a->ccwSucc;
    }
    else {
      if (a->poliDer != NULL &&
	  (d1 = ((x - a->poliDer->x) * (x - a->poliDer->x) + 
		 (y - a->poliDer->y) * (y - a->poliDer->y))) < d2) {
	d2 = d1;
	escogido = a->poliDer;
      }
      a = a->ccwPred;
    }
  } while(a != polIni->a);
  
  if (escogido == NULL) {
    return polIni;
  }
  return encontrarPoligono2(escogido, x, y);
}

void insertarSitio (MallaVoronoi *m, double *x1, double *y1, double valor) {
  PoligonoVoronoi *pol, *pol1, *pol2, *pol3;
  PuntoVoronoi *q, *qSel = NULL;
  AristaVoronoi *a;
  Lista *TPuntos, *TAristas, *aristasBorde;
  double H1, HMin = 1e300;
  NodoLista *n;
  double x, y;
  x = (*x1);
  y = (*y1);
  
  if (!validarMallaPoligonosAristas(m)) {
    fprintf (stderr, "ERROR al comienzo de insertarSitio\n");
    fprintf (stderr, "      insertando (%g, %g, %g)\n", x, y, valor);
    exit (1);
  }
  //printf ("Buscando poligono\n");
  pol1 = encontrarPoligono (m, x, y);
  if (DEBUG_MALLA) 
    printf ("Poligono encontrado: (%g, %g)\n", pol1->x, pol1->y);
  if (pol1 == NULL) {
    fprintf(stderr, 
	    "ERROR en insertarSitio: no se pudo encontrar el poligono.\n");
    return;
  }
  //if ((pol1->x == x) && (pol1->y == y)) {
  if (pow(pol1->x - x, 2) + pow(pol1->y - y, 2) < pow(PRECISION_DOUBLE, 2)) {
    //fprintf(stderr, "ERROR en insertarSitio: el Poligono (%g, %g) ya existe.\n",
    //x, y);
    (*x1) += (double) (random() / RAND_MAX - 0.5) * 1e-3;
    (*y1) += (double) (random() / RAND_MAX - 0.5) * 1e-3;
    fprintf(stderr, "      nuevo poligono = (%g, %g).\n", (*x1), (*y1));
    insertarSitio (m, x1, y1, valor);
    return;
  }
  
  if (DEBUG_MALLA) {
    printf ("-----insertando (%g, %g, %g)\n", x, y, valor);
  }
  pol = newPoligonoVoronoi(x, y, valor, m->nPols++);
  n = newNodoLista (pol, NULL, NULL);
  pol->n = n;
  insertarNodo (m->poligonos, n);
  
  // Buscamos el punto con el menor valor de H
  if ((a = pol1->a) == NULL) {
    fprintf(stderr, 
	    "ERROR en insertarSitio: no se pudo encontrar la arista inicial.\n");
    return;
  }
  
  if (a->poliDer == pol1) {
    a = a->cwPred;
  }
  else {
    a = a->cwSucc;
  }
  // printf("Arista posible error: \n");
  //imprimirArista(a);
  if (a->poliDer == pol1) {
    pol3 = a->poliIzq;
  }
  else if (a->poliIzq == pol1){
    pol3 = a->poliDer;
  }
  else {
    fprintf (stderr, "ERROR1 en insertarSitio: a no calza\n");
    fprintf (stderr, "       a->poliDer = (%g, %g)\n", 
	     a->poliDer->x, a->poliDer->y);
    fprintf (stderr, "       a->poliIzq = (%g, %g)\n", 
	     a->poliIzq->x, a->poliIzq->y);
    fprintf (stderr, "       pol1 = (%g, %g)\n", 
	     pol1->x, pol1->y);
    exit (1);
  }
  a = pol1->a;

  //printf (">>pol1 = (%g, %g)\n", pol1->x, pol1->y);

  do {
    pol2 = pol3;
    //printf (">>a: ");
    //imprimirArista (a);
    if (a->poliDer == pol1) {
      pol3 = a->poliIzq;
      q = a->ptoIni;
      a = a->ccwSucc;
    }
    else if (a->poliIzq == pol1) {
      pol3 = a->poliDer;
      q = a->ptoFin;
      a = a->ccwPred;
    }
    else {
      fprintf (stderr, "ERROR2 en insertarSitio: a no calza\n");
      fprintf (stderr, "       a->poliDer = (%g, %g)\n", 
	       a->poliDer->x, a->poliDer->y);
      fprintf (stderr, "       a->poliIzq = (%g, %g)\n", 
	       a->poliIzq->x, a->poliIzq->y);
      fprintf (stderr, "       pol1 = (%g, %g)\n", 
	       pol1->x, pol1->y);
      exit (1);
    }
    //printf (">>q = (%g, %g)\n", q->x, q->y);
    if (pol2 != NULL && pol3 != NULL && (H1 = H(pol1, pol3, pol2, pol)) < HMin) {
      //printf ("--pol1 = (%g, %g)\n", pol1->x, pol1->y);
      //printf ("--pol2 = (%g, %g)\n", pol2->x, pol2->y);
      //printf ("--pol3 = (%g, %g)\n", pol3->x, pol3->y);
      //printf ("-- H = %g\n", H1);
      //printf ("--q = (%g, %g)\n", q->x, q->y);
      HMin = H1;
      qSel = q;
    }
  } while (a != pol1->a);

  TPuntos = newLista ();
  //insertarNodo (TPuntos, newNodoLista(qSel, NULL, NULL));
  TAristas = newLista();
  aristasBorde = newLista();
  if (DEBUG_MALLA) {
    printf ("--qSel = (%g, %g)\n   a: ", qSel->x, qSel->y);
    if (qSel->a != NULL) {
      imprimirArista (qSel->a);
    }
    else {
      printf ("NULL\n");
    }
  }
  buscarPuntos (qSel, NULL, pol, TPuntos, TAristas, aristasBorde);
  if (DEBUG_MALLA) printf ("puntos encontrados\n");

  cambiarPoligonosAristasAEliminar (TAristas);

  if (aristasBorde->primero == NULL) {
    fprintf (stderr, 
	     "ERROR en insertarSitio (%g, %g, %g): aristasBorde->primero == NULL\n", *x1, *y1, valor);
    fprintf (stderr, "HMin = %g\n", HMin);
    imprimirMallaArchivo (m, "ERROR_NULL.dat");
  }
  //printf ("ok\n");
  //aIni = (AristaVoronoi *) aristasBorde->primero->info;
  //imprimirArista(aIni);

  crearVertices (m, aristasBorde, pol, pol1);

  if (!validarMallaPoligonosAristas(m)) {
    fprintf (stderr, "ERROR despues de crearVertices en insertarSitio\n");
    fprintf (stderr, "      insertando (%g, %g, %g)\n", x, y, valor);
    exit (1);
  }

  for (n = TAristas->primero; n != NULL; n = n->sgte) {
    a = (AristaVoronoi *) n->info;
    eliminarNodo (m->aristas, a->n);
    eliminarAristaVoronoi (a);
  }
  eliminarLista (TAristas);

  for (n = TPuntos->primero; n != NULL; n = n->sgte) {
    q = (PuntoVoronoi *) n->info;
    eliminarNodo (m->puntos, q->n);
    free (q);
  }
  eliminarLista (TPuntos);

  eliminarLista (aristasBorde);

  if (!validarMallaPoligonosAristas(m)) {
    fprintf (stderr, "ERROR al finalizar insertarSitio\n");
    fprintf (stderr, "      insertando (%g, %g, %g)\n", x, y, valor);
    exit (1);
  }
  //printf ("Insertado (%g, %g, %g)\n", pol->x, pol->y, pol->valor);
}

/* Cambia las aristas a las que apuntan los poligonos de las
   aristas a eliminar */
void cambiarPoligonosAristasAEliminar (Lista *aristas) {
  NodoLista *n;
  AristaVoronoi *a;

  for (n = aristas->primero; n != NULL; n = n->sgte) {
    a = (AristaVoronoi *) n->info;
    if (a->poliDer->a == a) {
      cambiarAristaPoligono (a->poliDer, aristas);
    }
    if (a->poliIzq->a == a) {
      cambiarAristaPoligono (a->poliIzq, aristas);
    }
  }
}

/* Cambia la arista a la que apunta pol por otra que no este en 
   la lista */
void cambiarAristaPoligono (PoligonoVoronoi *pol, 
			    Lista *aristas){
  AristaVoronoi *a;
  NodoLista *n;
  int seEncontro;
  
  a = pol->a;
  do {
    seEncontro = FALSE;
    for (n = aristas->primero; n != NULL && !seEncontro; n = n->sgte) {
      if (a == (AristaVoronoi *) n->info) {
	seEncontro = TRUE;
      }
    }
    if (!seEncontro) {
      pol->a = a;
      return;
    }
    if (a->poliDer == pol) {
      a = a->ccwSucc;
    }
    else {
      a = a->ccwPred;
    }
  } while(a != pol->a);
  fprintf (stderr, "ERROR en cambiarAristaPoligono: no se puede cambiar la arista de pol\n");
}

double H (PoligonoVoronoi *pi, PoligonoVoronoi *pj,
	  PoligonoVoronoi *pk, PoligonoVoronoi *p) {
  double J2, J3, J4, ik2, jk2, lk2;
  
  ik2 = pow(pi->x - pk->x, 2) + pow(pi->y - pk->y, 2);
  jk2 = pow(pj->x - pk->x, 2) + pow(pj->y - pk->y, 2);
  lk2 = pow(p->x - pk->x, 2) + pow(p->y - pk->y, 2);
  
  J2 = (pi->y - pk->y) * 0.5 * (jk2) - (pj->y - pk->y) * 0.5 * (ik2);
  J3 = (pi->x - pk->x) * 0.5 * (jk2) - (pj->x - pk->x) * 0.5 * (ik2);
  J4 = (pi->x - pk->x) * (pj->y - pk->y) - (pj->x - pk->x) * (pi->y - pk->y);
  
  return J2 * (p->x - pk->x) - J3 * (p->y - pk->y) + 0.5 * J4 * lk2;
}

/* Metodo que busca los puntos que caen dentro de pol comenzando en q.
 * Retorna TRUE en caso de que q caiga dentro de pol.
 */
int buscarPuntos (PuntoVoronoi *q, PuntoVoronoi *qPadre,
		  PoligonoVoronoi *pol, Lista *TPuntos, 
		  Lista *aristasDentro, Lista *aristasBorde) {

  AristaVoronoi *a;
  PuntoVoronoi *q1;
  PoligonoVoronoi *pol1, *pol2, *pol3;
  int i, dentro;

  if (DEBUG_MALLA) printf ("--- buscando desde punto (%g, %g), cuyo padre es ", q->x, q->y);
  if (DEBUG_MALLA) {
    if (qPadre != NULL) printf ("(%g, %g)\n", qPadre->x, qPadre->y);
    else printf ("NULL\n");
  }
  a = q->a;
  if (a->ptoIni == q) {
    pol1 = a->poliDer;
    pol2 = a->poliIzq;
  }
  else {
    pol2 = a->poliDer;
    pol1 = a->poliIzq;    
  }

  //printf("ok1\n");
  //imprimirArista (a);
  //printf("ok2\n");
  if (a->ptoFin == q) {
    q1 = a->ptoIni;
    pol3 = (a->ccwSucc->ptoIni == q)? a->ccwSucc->poliIzq : a->ccwSucc->poliDer;
  }
  else {
    q1 = a->ptoFin;
    pol3 = (a->ccwPred->ptoIni == q)? a->ccwPred->poliIzq : a->ccwPred->poliDer;
    //printf ("a->ccwPred->poliDer: (%g, %g)\n",
    //a->ccwPred->poliDer->x, a->ccwPred->poliDer->y);
  }
  
  if (pol1 == NULL || pol2 == NULL || pol3 == NULL) {
    return FALSE;
  }
  /*
  printf("poligonos para (%g, %g): \n", q->x, q->y);
  printf(" pol1 = (%g, %g)\n", pol1->x, pol1->y);
  printf(" pol2 = (%g, %g)\n", pol2->x, pol2->y);
  printf(" pol3 = (%g, %g)\n", pol3->x, pol3->y);
  printf(" pol  = (%g, %g)\n", pol->x, pol->y);
  printf("H : %g\n", H (pol1, pol2, pol3, pol));
  */

  if (H (pol1, pol2, pol3, pol) < -PRECISION_DOUBLE) {
    if (DEBUG_MALLA) 
      printf ("----insertando punto (%g, %g), H = %g\n", q->x, q->y,
	      H (pol1, pol2, pol3, pol));
    insertarNodo(TPuntos, newNodoLista(q, NULL, NULL));
    
    for (i = 0; i < 3; i++) {
      if (q1 != qPadre) {
	dentro = buscarPuntos(q1, q, pol, TPuntos, 
			      aristasDentro, aristasBorde);
	if (dentro) {
	  insertarNodo(aristasDentro, newNodoLista(a, NULL, NULL));
	}
	else {
	  insertarNodo(aristasBorde, newNodoLista(a, NULL, NULL));
	}
      }
      a = (a->ptoFin == q)? a->ccwSucc : a->ccwPred;
      q1 = (a->ptoFin == q)? a->ptoIni : a->ptoFin;
    }
    return TRUE;
  }
  if (DEBUG_MALLA) 
    printf ("----punto (%g, %g) no insertado H = %g\n", q->x, q->y,
	    H (pol1, pol2, pol3, pol));
  if (qPadre == NULL) {
    //double HError = H (pol1, pol2, pol3, pol);
    printf ("ERROR: no se encontro estructura dentro del poligono\n");
  }
  return FALSE;
}

/*  Metodo para crear los nuevos vertices.
 *  pol es el nuevo poligono, y pol1 el poligono mas cercano.
 */
void crearVertices (MallaVoronoi *m, Lista *aristasBorde, 
		    PoligonoVoronoi *pol, PoligonoVoronoi *pol1) {
  NodoLista *n1, *n;
  PuntoVoronoi *q, *q1;
  AristaVoronoi *a1, *a2, *aIni, *a;

  aIni = (AristaVoronoi *) aristasBorde->primero->info;
  //printf ("oki1\n");
  if (puntoDentroPoligono (aIni, pol) == NULL) {
    imprimirMallaArchivo (m, "ERROR_puntoDentroPoligono.dat");
    exit(1);
  }
  q1 = (puntoDentroPoligono (aIni, pol) == aIni->ptoIni)? 
    aIni->ptoFin : aIni->ptoIni;


  // Creamos los nuevos puntos.
  for (n1 = aristasBorde->primero; n1 != NULL; n1 = n1->sgte) {
    //printf ("oki2\n");
    a1 = (AristaVoronoi *) n1->info;
    q = puntoDentroPoligono(a1, pol);
    if (q == NULL) {
      imprimirMallaArchivo (m, "ERROR_puntoDentroPoligono.dat");
      exit(1);
    }

    if (a1->ptoIni == q) {
      a1->ptoIni = circuncentro (pol, a1->poliDer, a1->poliIzq);
      n = newNodoLista (a1->ptoIni, NULL, NULL);
      a1->ptoIni->n = n;
      insertarNodo (m->puntos, n);

      a1->ptoIni->a = a1;
    }
    else {
      a1->ptoFin = circuncentro (pol, a1->poliDer, a1->poliIzq);
      n = newNodoLista (a1->ptoFin, NULL, NULL);
      a1->ptoFin->n = n;
      insertarNodo (m->puntos, n);

      a1->ptoFin->a = a1;
    }
  }

  // Creamos las nuevas aristas.
  a1 = aIni;
  q = (q1 == a1->ptoIni)? a1->ptoFin : a1->ptoIni;
  //q = puntoDentroPoligono(a1, pol);
  /*  q1 = a1->ptoIni;
  //if (a1 == NULL)
  //printf("oki\n");
  if (pow (q->x - pol->x, 2) + pow (q->y - pol->y, 2) > 
      pow (q1->x - pol->x, 2) + pow (q1->y - pol->y, 2)) {
    q = q1;
    q1 = a1->ptoFin;
    }*/
  
  do {
    a2 = siguienteArista (a1, q, aristasBorde);
    if (a1->ptoFin == q) {
      if (a1->poliDer == a2->poliDer) {
	a = newAristaVoronoi (a1->ptoFin, a2->ptoIni);
	a2->cwPred = a;
	q = a2->ptoIni;
      }
      else {
	a = newAristaVoronoi (a1->ptoFin, a2->ptoFin);
	a2->cwSucc = a;
	q = a2->ptoFin;
      }
      a->poliDer = a1->poliDer;
      a1->ccwSucc = a;
    }
    else {
      if (a1->poliIzq == a2->poliIzq) {
	a = newAristaVoronoi (a1->ptoIni, a2->ptoFin);
	a2->cwSucc = a;
	q = a2->ptoFin;
      }
      else {
	a = newAristaVoronoi (a1->ptoIni, a2->ptoIni);
	a2->cwPred = a;
	q = a2->ptoIni;
      }
      a->poliDer = a1->poliIzq;      
      a1->ccwPred = a;
    }
    a->poliIzq = pol;
    a->cwPred = a1;
    a->ccwSucc = a2;
    calcularSinCos(a);
    n1 = newNodoLista(a, NULL, NULL);
    a->n = n1;
    insertarNodo(m->aristas, n1);
    //printf ("------------");
    //imprimirArista(a2);
    
    a1 = a2;

  } while (a1 != aIni);

  // Terminamos de conectar las nuevas aristas (interior a pol).
  //aIni = a1;
  //q = puntoDentroPoligono(a1, pol);
    /*a1->ptoFin;
  q1 = a1->ptoIni;
  if (pow (q->x - pol->x, 2) + pow (q->y - pol->y, 2) > 
      pow (q1->x - pol->x, 2) + pow (q1->y - pol->y, 2)) {
    q = q1;
    q1 = a1->ptoFin;
    }*/
  
  do {
    if (a1->ptoFin == q) {
      a2 = a1->ccwSucc;
      a = a1->cwSucc;
    }
    else {
      a2 = a1->ccwPred;
      a = a1->cwPred;
    }
    a2->ccwPred = a;
    a->cwSucc = a2;
    
    q = a2->ptoFin;
    a1 = a2->ccwSucc;
    /*
    if (a->ptoIni != a->cwPred->ptoIni && 
	a->ptoIni != a->cwPred->ptoFin) {
      fprintf (stderr, "ERROR en crearVertices: a->ptoIni no esta en cwPred.\n");
      printf ("     a : ");
      imprimirArista (a);
      printf ("     a->cwPred: ");
      imprimirArista (a->cwPred);	
      exit(1);
    }
    if (a->ptoIni != a->ccwPred->ptoIni && 
	a->ptoIni != a->ccwPred->ptoFin) {
      fprintf (stderr, "ERROR en crearVertices: a->ptoIni no esta en ccwPred.\n");
      printf ("     a : ");
      imprimirArista (a);
      printf ("     a->ccwPred: ");
      imprimirArista (a->ccwPred);	
      exit (1);
      }
    if (a->ptoFin != a->cwSucc->ptoIni && 
	a->ptoFin != a->cwSucc->ptoFin) {
      fprintf (stderr, "ERROR en crearVertices: a->ptoFin no esta en cwSucc.\n");
      printf ("     a : ");
      imprimirArista (a);
      printf ("     a->cwSuc: ");
      imprimirArista (a->cwSucc);	
      exit (1);
    }
    if (a->ptoFin != a->ccwSucc->ptoIni && 
	a->ptoFin != a->ccwSucc->ptoFin) {
      fprintf (stderr, "ERROR en crearVertices: a->ptoFin no esta en ccwSucc.\n");
      printf ("     a : ");
      imprimirArista (a);
      printf ("     a->ccwSucc: ");
      imprimirArista (a->ccwSucc);	
      exit (1);
    }
    */
  } while (a1 != aIni);
 
  pol->a = a;
}

/* Funcion que retorna el pto de a que cae dentro del poligono
   pol.
*/

PuntoVoronoi *puntoDentroPoligono (AristaVoronoi *a,
				   PoligonoVoronoi *pol) {
  //double dx, dy, y0, xc, yc;
  PoligonoVoronoi *pol1, *pol2, *pol3;
  //AristaVoronoi *a_aux;
  double HIni, HFin;
  
  if (DEBUG_MALLA) {
    printf ("buscando punto de ");
    imprimirArista (a);
    printf ("que caiga en (%g, %g)\n ", pol->x, pol->y);
  }
  
  // Probando a->ptoIni;
  pol1 = a->poliDer;
  pol2 = a->poliIzq;
  if (a->cwPred->poliDer == a->poliDer) {
    pol3 = a->cwPred->poliIzq;
  }
  else {
    pol3 = a->cwPred->poliDer;
  }
  
  if (DEBUG_MALLA) {
    if (pol1 != NULL) printf ("1) pol1 = (%g, %g)\n", pol1->x, pol1->y);
    else printf ("1) pol1 = NULL\n");
    if (pol2 != NULL) printf ("1) pol2 = (%g, %g)\n", pol2->x, pol2->y);
    else printf ("1) pol2 = NULL\n");
    if (pol3 != NULL) printf ("1) pol3 = (%g, %g)\n", pol3->x, pol3->y);
    else printf ("1) pol3 = NULL\n");
  }

  if (DEBUG_MALLA) 
    if (pol1 != NULL && pol2 != NULL && pol3 != NULL)
      printf ("H1 = %g\n", H (pol1, pol2, pol3, pol));
  
  if (pol1 != NULL && pol2 != NULL && pol3 != NULL) {
    HIni = H (pol1, pol2, pol3, pol); 
    //printf ("punto de ");
    //imprimirArista(a);
    //printf (" dentro de (%g, %g) = (%g, %g)\n", 
    //    pol->x, pol->y, a->ptoIni->x, a->ptoIni->y);
    //return a->ptoIni;
  }
  else {
    fprintf (stderr, 
	     "ERROR en mallaVoronoi::puntoDentroPoligono, pol1, pol2 o pol3 == NULL\n");
    exit (1);
  }

  // Probando a->ptoFin;
  pol1 = a->poliDer;
  pol3 = a->poliIzq;
  if (a->cwSucc->poliDer == a->poliIzq) {
    pol2 = a->cwSucc->poliIzq;
  }
  else {
    pol2 = a->cwSucc->poliDer;
  }
  if (DEBUG_MALLA) {
    if (pol1 != NULL) printf ("2) pol1 = (%g, %g)\n", pol1->x, pol1->y);
    else printf ("2) pol1 = NULL\n");
    if (pol2 != NULL) printf ("2) pol2 = (%g, %g)\n", pol2->x, pol2->y);
    else printf ("2) pol2 = NULL\n");
    if (pol3 != NULL) printf ("2) pol3 = (%g, %g)\n", pol3->x, pol3->y);
    else printf ("2) pol3 = NULL\n");
  }

  if (DEBUG_MALLA) 
    if (pol1 != NULL && pol2 != NULL && pol3 != NULL)
      printf ("H2 = %g\n", H (pol1, pol2, pol3, pol));
  
  if (pol1 != NULL && pol2 != NULL && pol3 != NULL) {
    HFin = H (pol1, pol2, pol3, pol);
    if (HFin < -PRECISION_DOUBLE && HFin < HIni) {
    //printf ("punto de ");
    //imprimirArista(a);
    //printf (" dentro de (%g, %g) = (%g, %g)\n", 
    //    pol->x, pol->y, a->ptoFin->x, a->ptoFin->y);
    return a->ptoFin;
    }
  }

  if (HIni < -PRECISION_DOUBLE) {
    return a->ptoIni;
  }

  fprintf (stderr, 
	   "ERROR puntoDentroPoligono: no se encuentra el punto\n");
  printf ( "      Arista ");
  imprimirArista(a);
  printf ( "      pol (%g, %g, %g)\n", pol->x, pol->y, pol->valor);

  return NULL;
  //exit (1);
  

  /* Version antes del 27/10/2005
  if (a->poliDer != NULL) {
    pol1 = a->poliDer;
  }
  else if (a->poliIzq != NULL) {
    pol1 = a->poliIzq;
  }
  else {
    fprintf (stderr, 
	     "ERROR puntoDentroPoligono: arista solo tiene poligonos NULL");
    exit (1);
  }

  // punto en el centro de pol y pol1
  xc = (pol1->x + pol->x) / 2;
  yc = (pol1->y + pol->y) / 2;
  
  dx = - (pol1->y - pol->y);
  dy = pol1->x - pol->x;
  //  if (dy == 0) {
  //dy = 1e-100;
    //return NULL;
  //}
  y0 = yc - dy * xc / dx;
  
  if (dy * (pol->x - xc) + dx * (yc - pol->y) < 0) {
    return (dy * (a->ptoIni->x - xc) + dx * (yc - a->ptoIni->y) < 0)? 
      a->ptoIni : a->ptoFin; 
  }
  else {
    return (dy * (a->ptoIni->x - xc) + dx * (yc - a->ptoIni->y) < 0)? 
      a->ptoFin : a->ptoIni;
  }
  */
}

/* Funcion que busca la siguiente arista a a1 perteneciente a aristasBorde
   en la direccion anti-horario. q es el punto de a1 dentro del nuevo 
   poligono.
*/
AristaVoronoi *siguienteArista(AristaVoronoi *a1, 
			       PuntoVoronoi *q,
			       Lista *aristasBorde) {
  AristaVoronoi *a2;
  NodoLista *n;

  if (q == a1->ptoIni) {
    a2 = a1->ccwPred;
  }
  else if (q == a1->ptoFin){
    a2 = a1->ccwSucc;
  }
  else {
    fprintf(stderr, "ERROR en siguienteArista: q no pertenece a a1.\n");
    return NULL;
  }
  
  for (n = aristasBorde->primero; n != NULL; n = n->sgte) {
    if (a2 == (AristaVoronoi *) n->info) {
      return a2;
    }
  }
  
  if (a2->cwPred == a1) {
    return siguienteArista (a2, a2->ptoFin, aristasBorde);
  }
  else if (a2->cwSucc == a1) {
    return siguienteArista (a2, a2->ptoIni, aristasBorde);
  }
  
  fprintf(stderr, "ERROR en siguienteArista: a2 no conecta con a1.\n");
  return NULL;
}

void imprimirAristasMalla (MallaVoronoi *m) {
  NodoLista *n;
  AristaVoronoi *a;
  int i = 0;
  
  printf ("imprimiendo aristas malla\n\n");
  for (n = m->aristas->primero; n != NULL; n = n->sgte) {
    i++;
    a = (AristaVoronoi *)n->info;
    printf("%d) ", i);
    imprimirArista (a);
    if (a->poliDer != NULL) {
      printf ("   poliDer: (%g, %g)\n", a->poliDer->x, a->poliDer->y);
    }
    else {
      printf ("   poliDer: NULL\n");
    }
    if (a->poliIzq != NULL) {
      printf ("   poliIzq: (%g, %g)\n", a->poliIzq->x, a->poliIzq->y);
    }
    else {
      printf ("   poliIzq: NULL\n");
    }
    printf ("   cwPred: ");
    imprimirArista (a->cwPred);
    printf ("   ccwPred: ");
    imprimirArista (a->ccwPred);
    printf ("   cwSucc: ");
    imprimirArista (a->cwSucc);
    printf ("   ccwSucc: ");
    imprimirArista (a->ccwSucc);
    printf ("\n");
  }
}

void imprimirAristasMallaArchivo (MallaVoronoi *m, char *nombre, 
				  int nx, int ny) {
  FILE *archivo = fopen (nombre, "w");
  NodoLista *n;
  AristaVoronoi *a;
  int nAristas = 0, i, j, **mask;
  double **im;
 
  mask = (int **) malloc(nx * sizeof(int *));
  for (i = 0; i < nx; i++) {
    mask[i] = (int *) malloc (ny * sizeof(int));
  }

  im = toImage (m, nx, ny, mask);

  fprintf (archivo, "%d\t%d\n", nx, ny);

  for (j = ny - 1; j >= 0; j--) {
    for (i = 0; i < nx; i++) {
      fprintf(archivo, "%g ", im[i][j]);
    }
    fprintf(archivo, "\n");
  }

  for (i = 0; i < nx; i++) {
    free (im[i]);
  }
  free (im);

  for (n = m->aristas->primero; n != NULL; n = n->sgte) {
    nAristas++;
  }
  
  fprintf (archivo, "%d\n", nAristas);
  for (n = m->aristas->primero; n != NULL; n = n->sgte) {
    a = (AristaVoronoi *) n->info;
    fprintf (archivo, "%g\t%g\t%g\t%g\n", 
	     a->ptoIni->x, a->ptoIni->y, a->ptoFin->x, a->ptoFin->y);
  }
  fclose (archivo);
}

void imprimirMallaArchivo (MallaVoronoi *m, char *nombre) {
  FILE *archivo = fopen (nombre, "w");
  NodoLista *n;
  PoligonoVoronoi *pol;
  int nPols = 0;
  
  // Contamos los poligonos
  for (n = m->poligonos->primero; n != NULL; n = n->sgte) {
    nPols++;
  }
  nPols -= 3;
  printf ("n=%d\n", nPols);
  fprintf (archivo, "%d\n", nPols * 3);

  // Imprimimos los poligonos
  for (n = m->poligonos->primero; n != NULL; n = n->sgte) {
    pol = (PoligonoVoronoi *) n->info;
    fprintf (archivo, "%g\t%g\t%g\n", 
	     pol->x, pol->y, pol->valor);
  }

  fclose (archivo);
}

int checkAristas (MallaVoronoi *m, double tol) {
  int ret = TRUE;
  double sinIzq, sinDer, cosDer, cosIzq;
  NodoLista *n;
  AristaVoronoi *a;

  for (n = m->aristas->primero; n != NULL; n = n->sgte) {
    a = (AristaVoronoi *) n->info;
    sinIzq = a->sinIzq;
    sinDer = a->sinDer;
    cosIzq = a->cosIzq;
    cosDer = a->cosDer;
    calcularSinCos (a);
    if ( fabs (sinIzq - a->sinIzq) > tol ||
	 fabs (sinDer - a->sinDer) > tol ||
	 fabs (cosIzq - a->cosIzq) > tol ||
	 fabs (cosDer - a->cosDer) > tol) {
      printf ("ERROR en checkArista: ");
      imprimirArista (a);
      printf ("sinDer = %g, sinIzq = %g, cosDer = %g, cosIzq = %g vs \n",
	      a->sinDer, a->sinIzq, a->cosDer, a->cosIzq);
      printf ("sinDer = %g, sinIzq = %g, cosDer = %g, cosIzq = %g\n",
	      sinDer, sinIzq, cosDer, cosIzq);
      printf ("dif:     %g        %g        %g        %g\n\n",
	      a->sinDer - sinDer, a->sinIzq - sinIzq, a->cosDer - cosDer,
	      a->cosIzq - cosIzq);
      ret = FALSE;
    }
    if (a->poliDer != NULL && a->poliIzq != NULL && 
	(fabs (a->sinDer + a->sinIzq) > tol ||
	 fabs (a->cosDer + a->cosIzq) > tol)) {
      printf ("ERROR en checkArista: ");
      imprimirArista (a);
      printf ("sinDer = %g, sinIzq = %g, cosDer = %g, cosIzq = %g vs \n",
	      a->sinDer, a->sinIzq, a->cosDer, a->cosIzq);
      ret = FALSE;
    }
  }
  return ret;
}
/*
int round (double x){
  return floor (x + 0.5);
}
*/
int validarMallaPoligonosAristas (MallaVoronoi *m) {
  NodoLista *n;
  AristaVoronoi *a;
  PoligonoVoronoi *pol;

  for (n = m->poligonos->primero; n != NULL; n = n->sgte) {
    pol = (PoligonoVoronoi *) n->info;

    a = pol->a;
    do {
      if (a->poliDer == pol) {
	a = a->cwPred;
      }
      else if (a->poliIzq == pol) {
	a = a->cwSucc;
      }
      else {
	fprintf (stderr, "ERROR en validacion: arista no coincide con poligono.\n");
	imprimirMallaArchivo (m, "MallaValidacion.dat");
	fprintf (stderr, "      pol = (%g, %g, %g)\n", 
		 pol->x, pol->y, pol->valor);
	imprimirArista (a);
	/*
	if (a->poliDer == NULL) {
	  printf ("poliDer = NULL\n");
	}
	else {
	  printf ("poliDer = (%g, %g, %g)\n", 
		  a->poliDer->x, a->poliDer->y, a->poliDer->valor);
	}
	if (a->poliIzq == NULL) {
	  printf ("poliIzq = NULL\n");
	}
	else {
	  printf ("poliIzq = (%g, %g, %g)\n", 
		  a->poliIzq->x, a->poliIzq->y, a->poliIzq->valor);
	}
	*/
	return FALSE;
      }
      if (a->ptoIni != a->cwPred->ptoIni && 
	  a->ptoIni != a->cwPred->ptoFin) {
	fprintf (stderr, "ERROR en validacion: a->ptoIni no esta en cwPred.\n");
	printf ("     a : ");
	imprimirArista (a);
	printf ("     a->cwPred: ");
	imprimirArista (a->cwPred);	
	return FALSE;
      }
      if (a->ptoIni != a->ccwPred->ptoIni && 
	  a->ptoIni != a->ccwPred->ptoFin) {
	fprintf (stderr, "ERROR en validacion: a->ptoIni no esta en ccwPred.\n");
	printf ("     a : ");
	imprimirArista (a);
	printf ("     a->ccwPred: ");
	imprimirArista (a->ccwPred);	
	return FALSE;
      }
      if (a->ptoFin != a->cwSucc->ptoIni && 
	  a->ptoFin != a->cwSucc->ptoFin) {
	fprintf (stderr, "ERROR en validacion: a->ptoFin no esta en cwSucc.\n");
	printf ("     a : ");
	imprimirArista (a);
	printf ("     a->cwSuc: ");
	imprimirArista (a->cwSucc);	
	return FALSE;
      }
      if (a->ptoFin != a->ccwSucc->ptoIni && 
	  a->ptoFin != a->ccwSucc->ptoFin) {
	fprintf (stderr, "ERROR en validacion: a->ptoFin no esta en ccwSucc.\n");
	printf ("     a : ");
	imprimirArista (a);
	printf ("     a->ccwSucc: ");
	imprimirArista (a->ccwSucc);	
	return FALSE;
      }

    } while (a != pol->a);
  }
  return TRUE;
}

void ajustarPoligonoCuadrado (double *x, double *y) {
  double dy = ((*y) - 0.5), dx = ((*x) - 0.5), x_ori = *x, y_ori = *y;

  if ((*x) < 0) {
    (*x) = 0;
    (*y) = (*x - 0.5) * dy / dx + 0.5;
  }
  if ((*x) > 1) {
    (*x) = 1;
    (*y) = (*x - 0.5) * dy / dx + 0.5;
  }
  if ((*y) < 0) {
    (*y) = 0;
    (*x) = (*y - 0.5) * dx / dy + 0.5;
  }
  if ((*y) > 1) {
    (*y) = 1;
    (*x) = (*y - 0.5) * dx / dy + 0.5;
  }

  if (*x != x_ori) {
    printf ("Cambiado x de %g a %g ", x_ori, *x);
  }
  if (*y != y_ori) {
    printf ("y cambiado y de %g a %g\n", y_ori, *y);
  }
}
