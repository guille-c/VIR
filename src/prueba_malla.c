#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arista.h"
#include "lista.h"
#include "mallaVoronoi.h"

main() {
  MallaVoronoi *m;
  struct nodoLista *n;
  float im[40][40], *xv, *yv, *iv;
  int nx = 128, ny = 128, nv = 15, i, j;
  PoligonoVoronoi *pol1, *pol2, *pol3, *pol4;
  PuntoVoronoi *p;
  AristaVoronoi *a;
  struct lista *l = newLista();
  double x, y;
  FILE *archivo = fopen("imagen.txt", "w");
  /*
  xv = malloc(nv * sizeof(double));
  yv = malloc(nv * sizeof(double));
  iv = malloc(nv * sizeof(double));
  
  srand((unsigned int) time(0));
  for (i = 1; i <= nv; i++) {
    xv[i - 1] = (double) rand() / RAND_MAX;
    yv[i - 1] = (double) rand() / RAND_MAX;
    iv[i - 1] = i;
    printf ("(%g, %g) = %g\n", xv[i-1], yv[i-1], iv[i-1]);
  }
  submkim (xv, yv, iv, nv, nx, ny, im);

  /*
  a = newAristaVoronoi (newPuntoVoronoi (0.0, 0.0), 
			newPuntoVoronoi (1.0, 1.0));
  n = newNodoLista (a, NULL, NULL);
  insertarNodo(l, n);
  a->n = l->primero;
  insertarNodo(l, newNodoLista (newAristaVoronoi (newPuntoVoronoi (1.0, 1.0), 
						  newPuntoVoronoi (2.0, 2.0)),
				NULL, NULL));
  imprimirArista(a->n->info);
  imprimirArista(l->primero->info);

	
  /*	pol1 = newPoligonoVoronoi(0.0, 0.0, 1);
	pol2 = newPoligonoVoronoi(0.1, 0.0, 1);
	pol3 = newPoligonoVoronoi(0.0, 0.1, 1);
	pol4 = newPoligonoVoronoi(0.11, 0.11, 1);

	printf("H = %g\n", H(pol3, pol1, pol2, pol4));
  */
    
  m = newMallaVoronoi();

  //insertarSitio (m, 0.5, 0.6 - 1e-8, 1);
  //insertarSitio (m, 0.4, 0.7, 2);
  //insertarSitio (m, -0.05, 0.25, 3);
  //imprimirAristasMallaArchivo (m, "aristas.dat", nx, ny);

  
  srand((unsigned int) time(0));
  for (i = 1; i <= 10; i++) {
    insertarSitio(m, (double) rand() / RAND_MAX,  
		  (double) rand() / RAND_MAX, i);
  }
  
  imprimirAristasMallaArchivo (m, "10aristas.dat", nx, ny);
  for (i = 1; i <= 90; i++) {
    insertarSitio(m, (double) rand() / RAND_MAX,  
		  (double) rand() / RAND_MAX, i);
  }
  imprimirAristasMallaArchivo (m, "100aristas.dat", nx, ny);
  for (i = 1; i < 400; i++) {
    insertarSitio(m, (double) rand() / RAND_MAX,  
		  (double) rand() / RAND_MAX, i);
  }
  imprimirAristasMallaArchivo (m, "500aristas.dat", nx, ny);
  /*
  for (n = m->aristas->primero; n != NULL; n = n->sgte) {
    imprimirArista((AristaVoronoi *)n->info);
  }
  /*
  pol1 = encontrarPoligono(m, 0.49, 0.0);
  printf("poligono encontrado: (%g, %g)\n", pol1->x, pol1->y);
  imprimirAristasMallaArchivo (m, "aristas0.dat", nx, ny);
  
  printf ("\n--------------POLIGONO 1------------------------------------\n");
  insertarSitio(m, 0.5, 0.0, 1);
  imprimirAristasMallaArchivo (m, "aristas1.dat", nx, ny);
  //imprimirAristasMalla(m);
  printf ("\n--------------POLIGONO 2------------------------------------\n");
  insertarSitio(m, 0.0, 0.5, 2);
  imprimirAristasMallaArchivo (m, "aristas2.dat", nx, ny);
  //imprimirAristasMalla(m);
  printf ("\n--------------POLIGONO 3------------------------------------\n");
  insertarSitio(m, 1.0, 0.5, 3);
  imprimirAristasMallaArchivo (m, "aristas3.dat", nx, ny);
  imprimirAristasMalla(m);
  

  printf ("\n--------------POLIGONO 4------------------------------------\n");
  insertarSitio(m, 0.5, 1.0, 4);
  imprimirAristasMallaArchivo (m, "aristas4.dat", nx, ny);
  printf ("\n--------------POLIGONO 5------------------------------------\n");
  insertarSitio(m, 0.5, 0.5, 5);
  imprimirAristasMallaArchivo (m, "aristas5.dat", nx, ny);
  printf ("\n--------------POLIGONO 6------------------------------------\n");
  insertarSitio(m, 0.5, 0.75, 6);
  imprimirAristasMallaArchivo (m, "aristas6.dat", nx, ny);
  printf ("\n--------------POLIGONO 7------------------------------------\n");
  insertarSitio(m, 0.75, 0.75, 7);
  imprimirAristasMallaArchivo (m, "aristas7.dat", nx, ny);
  printf ("\n--------------POLIGONO 8------------------------------------\n");
  insertarSitio(m, 0.0, 0.0, 8);
  imprimirAristasMallaArchivo (m, "aristas8.dat", nx, ny);
  
  pol1 = m->poligonos->primero->info;
  if (dentro(pol1, 0.2499, 0.0)) {
    printf ("Dentro de (%g, %g)\n", pol1->x, pol1->y);
  }
  else {
    printf ("Fuera\n");
  }
  /*
  im = toImage(m, nx, ny);

  for (j = ny - 1; j >= 0; j--) {
    for (i = 0; i < nx; i++) {
      printf("%g ", im[i][j]);
    }
    printf("\n");
  }
/*
  x = 0.0 * 1.0 / 0.0;
  printf("%g\n", x);
/*  for (j = ny - 1; j >= 0; j--) {
    for (i = 0; i < nx; i++) {
      fprintf(archivo, "%g ", im[i][j]);
    }
    fprintf(archivo, "\n");
  }
/*
  pol1 = newPoligonoVoronoi(0.5, 3 * sqrt(2)/4 + 1.0, 1);
  pol2 = newPoligonoVoronoi(-3 * sqrt(6)/4 + 0.5, -3 * sqrt(2)/4 + 0.5, 1);
  pol3 = newPoligonoVoronoi(3 * sqrt(6)/4 + 0.5, -3 * sqrt(2)/4 + 0.5, 1);
  //pol1 = newPoligonoVoronoi(0.5, 1.0);
  //pol2 = newPoligonoVoronoi(0.0, 0.0);
  //pol3 = newPoligonoVoronoi(1.0, 0.0);
  p = circuncentro(pol1, pol2, pol3);
  
  printf("(%g, %g), (%g, %g), (%g, %g)\n", 
	 pol1->x, pol1->y,
	 pol2->x, pol2->y,
	 pol3->x, pol3->y);
  printf("(%g, %g)\n", p->x, p->y);

  y = 0.0;
  x = - (pol2->y - pol3->y) / (pol2->x - pol3->x) * (y - p->y) + p->x;
  printf("(%g, %g)\n", x, y);
*/  

  //a = newAristaVoronoi(p1, p2);
  //imprimirArista(a);
  
  //insertarNodo(l, newNodoLista((void *) newPuntoVoronoi(1.2, 1.0), NULL, NULL));
  //printf("%g\n", ((PuntoVoronoi *) l->primero->info)->x);
  //imprimirListaInt(l);

}
