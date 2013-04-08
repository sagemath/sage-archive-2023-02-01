/******************************************************

Copyright 2004, 2010 Fabien de Montgolfier
fm@liafa.jussieu.fr

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
**********************************************************/

/********************************************************

    DECOMPOSITION MODULAIRE DE GRAPHES NON-ORIENTES

Cet algorithme construit l'arbre de decomposition modulaire
d'un graphe donne sous forme d'une matrice d'adjacence.
Il s'effectue en temps O(m log n) temps et $O(m)$ espace.
Il est la concatenation de deux algorithmes distincts.
Le premier realise une permutation factorisante des sommets du graphe
(pour cette notion, cf these de Christian Capelle)
grace a une technique d'affinage de partitions (cf Habib, Paul & Viennot)
Le second construit l'arbre de decomposition modulaire a partir de cette
permutation, cf mon memoire de DEA
Montpellier, decembre 2000
********************************************************/

//#include "dm_english.h"
#include <stdio.h>
#include <stdlib.h>

#define DEBUG 0                        /* si 0 aucune sortie graphique!! */

/* dm.h definit les constantes FEUILLE, MODULE, etc...
ainsi que les structures noeud et fils. Les autres
structures n'ont pas a etre connues par les programmes
exterieurs et sont donc definies ici. */


/* un sommet du graphe
(utilise dans la premiere partie seulement, ainsi que ce qui suit)*/
typedef struct Sommet {
  int place;/* numero du sommet dans la NOUVELLE numerotation */
  int nom;  /* numero du sommet dans l'ANCIENNE numerotation */
  /* On a donc sigma(nom)=place */
  struct Sadj *adj;
  struct SClasse *classe;
} sommet;

/* liste d'adjacence d'un sommet, DOUBLEMENT chainee*/
typedef struct Sadj {
  struct Sommet *pointe;
  struct Sadj *suiv;
  struct Sadj *prec;
} sadj;

/* classe de la partition courante,
   organisees en liste chainnee selon l'ordre de la partition */
typedef struct SClasse {
    int debut;
    int fin;
    struct Sommet *firstpivot;
    int inpivot;                /*indice de la classe dans le tableau pivot */
    int inmodule;                /* (resp module); -1 si non present */
    int whereXa;                /* lie le couple X/Xa: vaut
                                   0 si X n'est actuellement lie a aucun Xa
                                   -1 si Xa est a gauche
                                   +1 si Xa est a droite */
    struct SClasse *suiv;        /* forment une liste chainee */
    struct SClasse *prec;        /*...doublement */
} sclasse;

/* plein de parametres statiques que algo1() donne a Raffine() */
typedef struct Info {
  sclasse **pivot;
  int *ipivot;
  sclasse **module;
  int *imodule;
  int *numclasse;
  int *n;
} info;

/* clef a deux entrees utilisee pour le tri lineaire
   represente l'arrete ij */
typedef struct Clef2{
  int i; //sommet pointeur
  int nom; // nom du sommet pointe
  int place; //place du sommet pointe
} clef2;

/*************************************************************
utilitaires
*************************************************************/
void *fabmalloc(size_t s)
/* malloc sans erreur */
{
  void *p;
  p=malloc(s);
  if(p==NULL)
    {
      perror("Erreur de malloc!\n");
      exit(1);
    }
  return p;
}

int min(int a, int b)
{
  return (a<b) ? a : b;
}

int max(int a, int b)
{
    return (a > b) ? a : b;
}
/**************************************************************
Premiere passe de l'algorithme: il s'agit de trouver une
permutation factorisante du graphe. Nous utilisons les
techniques de raffinement de partition. Tout cela est
explique dans l'article de Habib, Viennot & Paul, dont je ne
fais ici que transcrire le travail.
****************************************************************/

void printS(sommet ** S, int n)
{
  /* imprimme S selon S et selon les classes */
  int i;
  sclasse *s;

  for (s = S[0]->classe; s != NULL; s = s->suiv) {
    printf("[ ");
    for (i = s->debut; i <= s->fin; i++)
      printf("%i ", 1 + S[i]->nom);
    printf("] ");
  }
  printf("\n");
}

sclasse *nouvclasse(sclasse * un, sclasse * deux)
{
  /* cree une nouvelle classe et l'insere entre un et deux;
     on suppose que si un et deux sont pas NULL alors
     FORCEMENT un=deux->suiv */

  sclasse *nouv;
  nouv = (sclasse *) fabmalloc(sizeof(sclasse));
  nouv->whereXa = 0;
  nouv->inpivot = -1;
  nouv->inmodule = -1;
  nouv->firstpivot = NULL;
  nouv->prec = un;
  if (un != NULL)                /* accroche pas en tete de chaine */
    un->suiv = nouv;
  nouv->suiv = deux;
  if (deux != NULL)                /* pas en queue de chaine */
    deux->prec = nouv;

    /* debut et fin ne sont PAS initialises! */
  return nouv;
}

void permute(sommet ** S, int a, int b)
{
  /* transpose les sommets a et b dans S */
  /* ne touche pas aux classes! */
  sommet *tmp;
  S[a]->place = b;
  S[b]->place = a;
  tmp = S[a];
  S[a] = S[b];
  S[b] = tmp;
}

void Raffiner(sommet ** S, sommet * p, sommet * centre, info * I)
{
  /* melange raffiner, pivotset, insertright et addpivot */
  sadj *a;                        /* parcours l'adjacence du pivot */
  sommet *x;                        /* sommet quiva changer de classe */
  sclasse *X, *Xa;                /* x in X; Xa nouv classe de x */
  sclasse *Z;
  sclasse **pivot;
  sclasse **module;
  int *ipivot, *imodule, *numclasse, n;

  if (DEBUG)
    printf("Raffinage avec le pivot %i\n", 1 + p->nom);
  pivot = I->pivot;
  module = I->module;
  ipivot = I->ipivot;
  imodule = I->imodule;
  numclasse = I->numclasse;
  n = *(I->n);

  for (a = p->adj; a != NULL; a = a->suiv) {
    x = a->pointe;
    X = x->classe;
    if (X == p->classe)
      continue;                /* on raffine pas la classe du pivot! */

    if (X->whereXa == 0) {
      /* c'est la premiere fois qu'on trouve un x
         appartenant a X lors de cet appel a raffiner */

      if ((centre->place < x->place && x->place < p->place)
          || (p->place < x->place && x->place < centre->place)) {
        /* insere a gauche */
        Xa = nouvclasse(X->prec, X);
        (*numclasse)++;
        permute(S, x->place, X->debut);
        X->debut++;
        X->whereXa = -1;
        Xa->whereXa = 1;        /* besoin dans le second tour */
      }
      else {                /* insere a droite */

        Xa = nouvclasse(X, X->suiv);
        (*numclasse)++;
        permute(S, x->place, X->fin);
        X->fin--;
        X->whereXa = 1;
        Xa->whereXa = -1;
      }
      x->classe = Xa;
      Xa->debut = x->place;
      Xa->fin = x->place;
    }
    else {
      if (X->whereXa == -1) {
        Xa = X->prec;
        permute(S, x->place, X->debut);
        X->debut++;
        Xa->fin++;
      }
      else {
        Xa = X->suiv;
        permute(S, x->place, X->fin);
        X->fin--;
        Xa->debut--;
      }
      x->classe = Xa;
    }
  }

  for (a = p->adj; a != NULL; a = a->suiv)
    /* deuxieme couche! Maintenant on va faire les addpivot,
           et remettre les whereXa a 0
           Noter qu'on lit les Xa et plus les X */
    {
      x = a->pointe;
      Xa = x->classe;
      if (Xa->whereXa == 0)
        continue;                /* deja remis a zero! */
      if (Xa->whereXa == -1)
        X = Xa->prec;
      else
        X = Xa->suiv;

      if (X->debut > X->fin) {
        /*on a trop enleve! X est vide
          -> on va le supprimer mechamment */

        (*numclasse)--;
        if (Xa->whereXa == 1) {        /*deconnecte */
          Xa->suiv = X->suiv;
          if (Xa->suiv != NULL)
            Xa->suiv->prec = Xa;
        }
        else {
          Xa->prec = X->prec;
          if (Xa->prec != NULL)
            Xa->prec->suiv = Xa;
        }
        Xa->inpivot = X->inpivot;
        if (X->inpivot != -1)        /* ecrase X dans pivot */
          pivot[X->inpivot] = Xa;
        Xa->inmodule = X->inmodule;
        if (X->inmodule != -1)        /* ecrase X dans pivot */
          module[X->inmodule] = Xa;

        Xa->whereXa = 0;
        continue;
      }

      /* Maintenant on fait addpivot(X,Xa)
         noter que X et Xa sont non vides */

      if (X->inpivot == -1) {
        if ((X->inmodule != -1)
            && (X->fin - X->debut < Xa->fin - Xa->debut)) {
          /* remplace X par Xa dans module */
          module[X->inmodule] = Xa;
          Xa->inmodule = X->inmodule;
          X->inmodule = -1;
          if (DEBUG)
            printf("Dans module %i-%i ecrase %i-%i\n",
                   1 + S[Xa->debut]->nom, 1 + S[Xa->fin]->nom,
                   1 + S[X->debut]->nom, 1 + S[X->fin]->nom);
        }
        else {
          if (X->inmodule == -1) {
            if (X->fin - X->debut < Xa->fin - Xa->debut)
              Z = Xa;
            else
              Z = X;
            /* ajoute Z (=max(X,Xa)) a module */
            module[(*imodule)] = Z;
            Z->inmodule = (*imodule);
            (*imodule)++;
            if (DEBUG)
              printf("module empile:%i-%i\n",
                     1 + S[Z->debut]->nom, 1 + S[Z->fin]->nom);
          }
        }
      }

      if (X->inpivot != -1)
        Z = Xa;
      else if (X->fin - X->debut < Xa->fin - Xa->debut)
        Z = X;
      else
        Z = Xa;
      /* on empile Z dans pivot */
      pivot[(*ipivot)] = Z;
      Z->inpivot = (*ipivot);
      (*ipivot)++;
      if (DEBUG)
        printf("pivot empile: %i-%i\n", 1 + S[Z->debut]->nom,
               1 + S[Z->fin]->nom);
      X->whereXa = 0;
      Xa->whereXa = 0;
    }
  if (DEBUG) {
    printS(S, n);
    printf("\n");
  }
}

sommet **algo1(graphe G)
     /* Entree: un graphe G
        Sortie: une permutation factorisante de G,
        donnee sous la forme d'un tableau de structures Sommet ordonnees selon sigma.
        d'apres le travail de Habib/Paul/Viennot */
{
  int n; // nombre de sommets de G

  sclasse **pivot;                /*pile des pivots */
  int ipivot = 0;                /*indice sur la precedante */

    sclasse **module;                /*idem, modules */
    int imodule = 0;

    sclasse *singclasse;
    /*invariant: toute classe avant singclasse dans la chaine */
    /*a un seul element */
    int numclasse;                /* quand vaut n, on a fini!! */

    sclasse *C1;                /*premiere classe, tete de la chaine */
    sclasse *Y;                        /*classe qui raffine */
    sclasse *X;                        /*classe raffinee */
    sclasse *Xa, *Xc;                /* morceaux de X */
    sommet *x;                        /* x in X */
    sommet *y;                        /* y in Y */
    sommet *centre;                /* le centre du raffinage actuel */

    sommet **S;                        /*la permutation factorisante ! */

    int i, j;                        /*divers indices */
    sommet *scourant;                /* pour l'init */
    sadj *nextadj;                /*sommet adjacent suivant */
    adj *nextadj2;              /* idem mais de type adj */
    info Inf;                        /* diverses info a passer a raffiner */

    /* debut des initialisations */
    n=G.n;
    /*initialisation des tableaux */
    module = (sclasse **) fabmalloc(n * sizeof(sclasse *));
    pivot = (sclasse **) fabmalloc(n * sizeof(sclasse *));
    S = (sommet **) fabmalloc(n * sizeof(sommet *));
    /* on va initialiser la permutation factorisante,
       ainsi que chaque structure sommet */
    C1 = nouvclasse(NULL, NULL);
    numclasse = 1;
    singclasse = C1;
    C1->debut = 0;
    C1->fin = n - 1;
    for (i = 0; i < n; i++) {
        /* initialisation des sommets */
        /* notre bebe est le sommet i dans M */
        scourant = (sommet *) fabmalloc(sizeof(struct Sommet));
        scourant->nom = i;
        scourant->place = i;        /* a ce point S=identite */
        scourant->adj = NULL; /* pas encore d'adjacence */
        scourant->classe = C1;
        S[i] = scourant;
    }
    for (i = 0; i < n; i++)
      {
        nextadj2 = G.G[i];
        while(nextadj2 != NULL)
          {
            j=nextadj2->s; //numero du sommet pointe
            if((j<0)||(j>=n))
              {
                perror("Graphe invalide (numero de sommet erronne)!\n");
                exit(1);
              }
            nextadj = (sadj *) fabmalloc(sizeof(struct Sadj));
            //un nouveau sadj
            nextadj->pointe = S[j];
            nextadj->suiv = S[i]->adj; //tete de liste
            if(nextadj->suiv!=NULL)
              nextadj->suiv->prec=nextadj;
            nextadj->prec=NULL;
            S[i]->adj = nextadj;        /*et le tour est joue */
            nextadj2=nextadj2->suiv;
          }
        }
    /* NB: module et pivot sont vides */
    Inf.pivot = pivot;
    Inf.ipivot = &ipivot;
    Inf.module = module;
    Inf.imodule = &imodule;
    Inf.numclasse = &numclasse;
    Inf.n = &n;
    /* init terminnee */

    while (1) {
        while (ipivot > 0 || imodule > 0) {
            while (ipivot > 0) {
                /*cette boucle raffine selon tous les sommets
                   de la premiere classe dans pivot */

                Y = pivot[ipivot - 1];
                ipivot--;
                Y->inpivot = -1;

                for (i = Y->debut; i <= Y->fin; i++)
                    Raffiner(S, S[i], centre, &Inf);

                /* une optimisation de la fin de l'algo */
                if (numclasse == n)
                    return (S);
            }
            /*maintenant pivot est vide, mais peut-etre pas module */
            if (imodule > 0) {
                /* relance par un sommet (pas au pif...) */
                /* de chaque module qui le represente */
                Y = module[imodule - 1];
                imodule--;
                Y->inmodule = -1;
                y = S[Y->debut];        /* le firstpivot sera toujours... */
                Y->firstpivot = y;        /* le premier!! */
                if (DEBUG)
                    printf("module-pivot %i-%i: sommet %i\n",
                           1 + S[Y->debut]->nom, 1 + S[Y->fin]->nom,
                           1 + y->nom);
                Raffiner(S, y, centre, &Inf);
            }
        }
        /* a ce point, pivot et module sont vides...
           pas de pb! On va faire initpartition HERE */
        if (DEBUG)
            printf("\nInit Partition\n");
        /**** ajoute ici pour debbugger, mais moche!! */
        singclasse = S[0]->classe;
        while ((singclasse != NULL) &&
               (singclasse->debut == singclasse->fin))
          {
            singclasse = singclasse->suiv;
          }
        /* singclasse est la premiere classe
           non singlette, sauf si: */
        if (singclasse == NULL)
            /* on a n classes singlettes? ben c'est gagne! */
          {
            return (S);
          }
        if (singclasse == NULL && numclasse < n) {
            perror("c'est pas normal! Ca termine trop vite!\n");
            exit(1);
        }

        X = singclasse;
        x = X->firstpivot;
        if (x == NULL)
            x = S[X->debut];
        else                        /* remet firstpivot a NULL!! */
            X->firstpivot = NULL;

        if (DEBUG)
            printf("Relance dans le module %i-%i avec le sommet %i\n",
                   1 + S[X->debut]->nom, 1 + S[X->fin]->nom, 1 + x->nom);

        centre = x;                /*important! */
        /* astuce: on place {x} en tete de X
           ensuite, on raffine S selon x -> seule X est coupee
           il y a alors {x} X Xa
           -> on met {x} en queue de X et c'est bon!
           ainsi on a bien nonvoisins-x-voisons */
        Xc = nouvclasse(X->prec, X);
        numclasse++;
        x->classe = Xc;
        permute(S, x->place, X->debut);
        X->debut++;
        Xc->debut = x->place;
        Xc->fin = x->place;
        Raffiner(S, x, x, &Inf);
        /* X existe-il encore? */
        if (X->debut > X->fin)
            continue;
        /* echange de x et {x}. Init: -{x}-X- */
        Xc->suiv = X->suiv;
        if (X->suiv != NULL)
            X->suiv->prec = Xc;
        X->prec = Xc->prec;
        if (Xc->prec != NULL)
            Xc->prec->suiv = X;
        X->suiv = Xc;
        Xc->prec = X;
        permute(S, x->place, X->fin);
        Xc->debut = x->place;
        Xc->fin = x->place;
        X->debut--;
        X->fin--;
        //antibug?
        singclasse=X;
        /* now -X-{x}- */
        if (DEBUG)
            printS(S, n);
    }
}

/***************************************************************
Etape intermediaire: trier toutes les listes d'adjacence
selon S. ce sont les listes de type sadj qui sont concernees
***************************************************************/
int Calculm(graphe G)
/* compte le nombre d'arretes du graphe */
{
  int i,r; adj *a;
  r=0;
  for(i=0;i<G.n;i++)
    {
      a=G.G[i];
      while(a!=NULL)
        {
          a=a->suiv;
          r++;
        }
    }
  if(r%2!=0)
    {
      perror("Erreur: nombre impaire d'arrete, graphe non-oriente??\n");
      exit(1);
    }
  return r/2; // G symetrique!
}

void TrierTous(sommet **S, int n, int m)
/* trie chaque liste d'adjacence de S*/
{
  //n sommets, m arretes
  int i; // numero du sommet courant
  sadj *a,*atmp;// parcours sa liste d'adjacence
  clef2 *c; // enregistrement a trier
  int *tab1; clef2 **tab2; //tableaux du tri par seaux
  tab1=(int *)fabmalloc(n*sizeof(int));
  tab2=(clef2 **)fabmalloc(m * 2 * sizeof(clef2 *));
  for(i=0; i<n; i++)
    tab1[i]=0;

  // premiere passe: construit tab1:
  // tab[i] est la frequence du ieme (selon S) sommet
  for(i=0; i<n; i++)
    {
      a=S[i]->adj;
      while(a!=NULL)
        {
          tab1[i]++;
          a=a->suiv;
        }
    }
  //deuxieme passe: frequences cumulees a rebours
  // (car les listes d'adjacences se construisent a l'envers
  //tab1[n-1]--; // a cause des indices de tableau qui commence a zero
  //for(i=n-1;i>0;i--)
  //  tab1[i-1]+=tab1[i];

  //deuxieme passe: frequences cumulees
  for(i=1;i<n;i++)
    tab1[i]+=tab1[i-1];

  //troisieme passe: liste double
  for(i=0; i<n; i++)
    {
      a=S[i]->adj;
      while(a!=NULL)
        {
          /* cree un nouveau record */
          c=(clef2 *)fabmalloc(sizeof(struct Clef2));
          c->i=i;
          c->nom=a->pointe->nom;
          c->place=a->pointe->place;
          /* le place bien dans tab2 */
          tab1[c->place]--;
          tab2[tab1[c->place]]=c;
          /*et on continue */
          a=a->suiv;
        }
    }

  //quatrieme passe: detruit les vielles listes d'adjacence
  for(i=0; i<n; i++)
    {
      a=S[i]->adj;
      while(a!=NULL)
        {
          atmp=a->suiv;
          free(a);
          a=atmp;
        }
      S[i]->adj=NULL;
    }

  //derniere passe: reconstruit les listes d'adjacence
  for(i=0;i<2*m;i++)
    {
      c=tab2[i];
      a=(sadj *)fabmalloc(sizeof(struct Sadj));
      a->pointe=S[c->i];
      a->suiv=S[c->place]->adj; //insere en tete
      if(a->suiv!=NULL)
        a->suiv->prec=a;
      a->prec=NULL;
      S[c->place]->adj=a;
      //nettoie
      free(c);
   }
  free(tab1);
  free(tab2);
}


/***************************************************************
 Maintenant, la deuxieme partie de l'aglorithme
 On va, etant donne la matrice M construite a l'etape precedante,
 etablir l'arbre de decomposition modulaire.
 Tous les details sont dans mon memoire de DEA
****************************************************************/
noeud *nouvnoeud(int type, noeud * pere, int sommet, int n)
{
    /* cree un nouveau noeud. Noter que l'on est oblige de passer n
       comme parametre car les bords et separateurs droits doivent
       etre initilises avec des valeurs >n */
    noeud *nn;
    static int compteur = 0;
    /*pour donner un ID unique aux noeuds. juste pour debug */

    nn = (noeud *) fabmalloc(sizeof(noeud));
    nn->type = type;
    nn->pere = pere;
    /* nn->fpere ne peut etre deja mis a jour... */
    nn->sommet = sommet;
    nn->ps = n + 2;
    nn->ds = -2;
    /*ces valeurs pour distinguer "non calcule" (-2) */
    /*   de "abscence de separateur" (-1). De plus, on fera des min et des */
    /*   max sur les bords */
    nn->bg = n + 2;
    nn->bd = -2;                /* idem */

    nn->fils = NULL;
    nn->lastfils = NULL;
    nn->id = compteur;
    compteur++;
    return nn;
}

void ajoutfils(noeud * pere, noeud * nfils)
{
    fils *nf;
    /* noter que c'est un ajout en queue! */
    nf = (fils *) fabmalloc(sizeof(fils));
    nf->pointe = nfils;
    nf->suiv = NULL;
    if (pere->fils == NULL)
        pere->fils = nf;        /* on cree le premier fils */
    else
        pere->lastfils->suiv = nf;        /* on ajoute nf a la chaine */
    pere->lastfils = nf;
    nfils->pere = pere;                /* normalement: redondant,mais... */
    nfils->fpere = nf;
}

void fusionne(noeud * pere, noeud * artefact)
{
    /*fusionne un artefact a son pere.
       utilise le champ fpere qui permet de savoir ou se greffer
       une structure fils sera detruite dans l'operation: artefact->fils */
    fils *greffe;
    fils *f;
    /* met a jour la liste des peres */
    f = artefact->fils;
    while (f != NULL) {
        f->pointe->pere = pere;        /*avant c'etait ancien... */
        /* f->pointe->fpere est inchange */
        f = f->suiv;
    }
    /* greffe la liste */
    greffe = artefact->fpere;
    artefact->lastfils->suiv = greffe->suiv;
    greffe->pointe = artefact->fils->pointe;
    greffe->suiv = artefact->fils->suiv;
    artefact->fils->pointe->fpere = greffe;        /*artefact->fils a disparu */
    if (pere->lastfils == greffe)
        pere->lastfils = artefact->lastfils;
}

void
extraire(noeud * ancien, noeud * nouveau, fils * premier, fils * dernier)
{
    /* extrait la liste [premier...dernier] des fils de l'ancien noeud,
       et en fait la liste des fils du nouveau noeud */
    fils *nf;                        /* il faut une structure fils de plus */
    fils *f;                        /*indice de mise a jour */
    nf = (fils *) fabmalloc(sizeof(fils));
    nf->pointe = premier->pointe;
    nf->suiv = premier->suiv;
    premier->pointe->fpere = nf;
    nouveau->pere = ancien;
    nouveau->fils = nf;
    nouveau->lastfils = dernier;
    nouveau->bg = premier->pointe->bg;
    nouveau->bd = dernier->pointe->bd;
    nouveau->ps = premier->pointe->bg;        /* nouveau est suppose etre un */
    nouveau->ds = dernier->pointe->bd;        /* module, donc bords=separateurs! */
    if (ancien->lastfils == dernier)
        ancien->lastfils = premier;
    /* ecrase l'ancier premier */
    nouveau->fpere = premier;
    premier->pointe = nouveau;
    premier->suiv = dernier->suiv;
    /* met a jour dernier */
    dernier->suiv = NULL;
    /* met a jour la liste des peres */
    f = nf;
    while (f != dernier->suiv) {
        f->pointe->pere = nouveau;        /*avant c'etait ancien... */
        f->pointe->fpere = premier;
        f = f->suiv;
    }
}

void printnoeud(noeud * N, int level)
{
    /* imprime recursivement l'arbre par parcours en profondeur */
    fils *ffils;
    noeud *nfils;
    int i;
    ffils = N->fils;

    for (i = 0; i < level - 1; i++)
        printf("  |");
    if (N->pere == NULL)
        printf(" ");
    else
        printf("  +-");
    switch (N->type) {
    case UNKN:
        printf("Noeud\n");
        break;
    case MODULE:
        printf("Module\n");
        break;
    case ARTEFACT:
        printf("Artefact\n");
        break;
    case SERIE:
        printf("Serie \n");
        break;
    case PARALLELE:
        printf("Parallele \n");
        break;
    case PREMIER:
        printf("Premier \n");
        break;
    }

    do {
        nfils = ffils->pointe;
        if (nfils->type == FEUILLE) {
            for (i = 0; i < level; i++)
                printf("  |");
            printf("  +--");
            printf("%i\n", 1 + nfils->nom);
        }
        else {
            printnoeud(nfils, level + 1);
        }
        ffils = ffils->suiv;
    }
    while (ffils != NULL);
}

void printarbre(noeud * N)
{
    printnoeud(N, 0);
}

noeud *algo2(graphe G, sommet **S)
{
/* algorithme de decomposition modulaire, deuxieme passe
entree: le graphe G, et sa permutation factorisante S.
sortie: un pointeur sur un arbre de decomposition modulaire
*/
    /* debug: S n'est utilise que pour mettre vrainom a jour */
    int n; //nombre de sommets du graphe
    int *ouvrantes;                /* tableau du nombre de parentheses ouvrantes */
    /* ouvrante[i]=3 ssi i-1(((i  */
    /* ouvrante [0]=3: (((0 */

    int *fermantes;                /* idem fermantes[i]=2 ssi i)))i+1
                                   fermante [n-1]=2 ssi n))) */
    int *ps;                        /* ps[i]=premier separateur de (i,i+1) */
    int *ds;

    int i, j;                        /*indices de paires ou de sommets */

    sadj *a1, *a2;                /* parcours de liste d'adjacence */

    noeud *racine;                /*racine du pseudocoardre */
    noeud *courant, *nouveau;        /* noeud courant du pseudocoarbre */
    noeud **pileinterne;        /* pile des modules pour les passes 3,5,5 */
    int indicepileinterne = 0;        /*pointeur dans cette pile */
    int taillepileinterne;        /* taille de la pile apres la 2eme passe */

    int *adjii;                 /* adjii[i]=1 ssi S[i] et S[i+1] sont */
                                   /*                             adjacents */
    /*PROPHASE : initialisations */
    n=G.n;
    ouvrantes = (int *) fabmalloc(n * sizeof(int));
    fermantes = (int *) fabmalloc(n * sizeof(int));
    ps = (int *) fabmalloc(n * sizeof(int));
    ds = (int *) fabmalloc(n * sizeof(int));
    pileinterne = (noeud **) fabmalloc((2 * n + 4) * sizeof(noeud *));
    adjii= (int *) fabmalloc(n*sizeof(int));
    /*pas plus de 2n+4 noeuds internes dans le pseudocoarbre */
    for (i = 0; i < n; i++) {
      ouvrantes[i] = 0;
      fermantes[i] = 0;
      adjii[i]=0;
    }

    /* remplit adjii qui dit si S[i] adjacent a S[i+1] */
    for(i=0; i<n-1; i++)
      {
         a1=S[i]->adj;
         while((a1!=NULL)&&(a1->pointe->place != i+1))
           a1=a1->suiv;
         if( a1 == NULL)
           adjii[i]=0;
         else // a1->pointe->place==i+1, donc i adj i+1
           adjii[i]=1;
      }
    adjii[n-1]=0; //perfectionnisme

    /* PREMIERE PASSE
       on va parentheser la permutation factorisante.
       tout bonnement, on lit S et on cherche les separateurs;
       apres quoi ouvrantes et fermantes sont ecrites
       complexite: O(n^2) */

    ouvrantes[0] = 1;
    fermantes[n - 1] = 1;        /* parentheses des bords */

    for (i = 0; i < n - 1; i++) {
      /*recherche de ps(i,i+1) */
      a1=S[i]->adj;
      a2=S[i+1]->adj;
      while((a1!=NULL) && (a2!=NULL) && (a1->pointe->place<i) &&
            (a2->pointe->place<i) && (a1->pointe->place == a2->pointe->place))
        {
          a1=a1->suiv;
          a2=a2->suiv;
        }

      //arbre de decision complique pour trouver le premier separateur!
      if( ((a1==NULL) && (a2==NULL))
          ||((a1==NULL) &&(a2->pointe->place >= i))
          ||((a2==NULL) && (a1->pointe->place >= i))
          ||((a1!=NULL) && (a2!=NULL) && (a1->pointe->place >= i) && (a2->pointe->place >= i)))
        //pas de separateur
        ps[i]=i+1;
      else
        {
          if((a1==NULL) || (a1->pointe->place >= i))
            ps[i]=a2->pointe->place;
          else if((a2==NULL) || (a2->pointe->place >= i))
            ps[i]=a1->pointe->place;
          else
            {
              if((a1->suiv!=NULL)&&(a1->suiv->pointe->place == a2->pointe->place))
                ps[i]=a1->pointe->place;
              else if ((a2->suiv!=NULL)&&(a2->suiv->pointe->place == a1->pointe->place))
                ps[i]=a2->pointe->place;
              else
                ps[i]=min(a1->pointe->place , a2->pointe->place);
            }
          ouvrantes[ps[i]]++;        /* marque la fracture gauche, if any */
          fermantes[i]++;
        }
      if (DEBUG)
        printf("ps(%i,%i)=%i\n", i , i+1, ps[i]);

        /*recherche de ds(i,i+1)
          plus penible encore!*/
      a1=S[i]->adj;
      if(a1!=NULL) // se place en queue de liste.
        while(a1->suiv!=NULL)
          a1=a1->suiv;
      a2=S[i+1]->adj;
      if(a2!=NULL)
        while(a2->suiv!=NULL)
          a2=a2->suiv;
      while((a1!=NULL) && (a2!=NULL) && (a1->pointe->place > i+1) &&
            (a2->pointe->place > i+1) && (a1->pointe->place == a2->pointe->place))
        {
          a1=a1->prec;
          a2=a2->prec;
        }
      if( ((a1==NULL) && (a2==NULL))
          ||((a1==NULL) && (a2->pointe->place <= i+1))
          ||((a2==NULL) && (a1->pointe->place <= i+1))
          ||((a1!=NULL) && (a2!=NULL) && (a1->pointe->place <= i+1) && (a2->pointe->place <= i+1)))
        //pas de separateur
        ds[i]=i+1;
      else
        {
          if((a1==NULL) || (a1->pointe->place <= i+1))
            ds[i]=a2->pointe->place;
          else if((a2==NULL) || (a2->pointe->place <= i+1))
            ds[i]=a1->pointe->place;
          else
            {
              if((a1->prec!=NULL)&&(a1->prec->pointe->place == a2->pointe->place))
                ds[i]=a1->pointe->place;
              else if((a2->prec!=NULL)&&(a2->prec->pointe->place == a1->pointe->place))
                ds[i]=a2->pointe->place;
              else
                ds[i]=max(a1->pointe->place , a2->pointe->place);
            }


          //ds[i] = j;
        ouvrantes[i + 1]++;        /* marque la fracture gauche, if any */
        fermantes[ds[i]]++;        /* attention aux decalages d'indices */
        }
      if (DEBUG)
        printf("ds(%i,%i)=%i\n", i,i+1,ds[i]);
      //S[i]->nom + 1,               S[i + 1]->nom + 1, S[ds[i]]->nom + 1);
    }

    /*DEUXIEME PASSE: construction du pseudocoarbre */

    racine = nouvnoeud(UNKN, NULL, -1, n);
    courant = racine;
    for (i = 0; i < n; i++) {
        /*1: on lit des parentheses ouvrantes: descentes */
        for (j = 0; j < ouvrantes[i]; j++) {
            /*Descente vers un nouveau noeud */
            nouveau = nouvnoeud(UNKN, courant, -1, n);
            ajoutfils(courant, nouveau);        /*on l'ajoute... */
            courant = nouveau;        /* et on descent */
            if (DEBUG)
                printf("(");
        }
        /* 2: on lit le sommet: feuille */
        nouveau = nouvnoeud(FEUILLE, courant, i, n);
        ajoutfils(courant, nouveau);        /*on ajoute le bebe... */
        (*nouveau).ps = i;
        (*nouveau).ds = i;
        (*nouveau).bg = i;
        (*nouveau).bd = i;
        nouveau->nom = S[i]->nom;
        /* et pourquoi pas i ? Je m'embrouille... */
        if (DEBUG)
            printf(" %i ", S[i]->nom + 1);

        /*3: on lit des parentheses fermantes: remontees */
        for (j = 0; j < fermantes[i]; j++) {
            /*ASTUCE: ici on va en profiter pour supprimer les
               noeuds a un fils, afin d'economiser une passe */
            if (courant->fils == courant->lastfils) {        /*just one son */
                courant->pere->lastfils->pointe = courant->fils->pointe;
                courant->fils->pointe->pere = courant->pere;
                courant->fils->pointe->fpere = courant->pere->lastfils;
                /* explication: le dernier fils de courant.pere est
                   actuellement courant himself. Il est donc pointe par
                   courant.pere.lastfils.pointe. Il suffit de changer ce
                   pointeur pour qu'il pointe maintenant non plus sur courant,
                   mais sur l'unique fils de courant: courant.fils.pointe.
                   Ouf! */
                /* NB: courant est maintenant deconnecte de l'arbre.
                   on pourrait faire un free() mais bon... */
            }
            else {
                /*on empile ce noeud interne.
                   L'ordre est celui de la postvisite d'un DFS */
                pileinterne[indicepileinterne] = courant;
                indicepileinterne++;
            }
            /* et dans tous les cas: on remonte! */
            courant = courant->pere;
            if (DEBUG)
                printf(")");
        }
    }
    if (DEBUG)
        printf("\n");

    /* on enleve un ptit defaut */
    racine = racine->fils->pointe;
    racine->pere = NULL;
    racine->fpere = NULL;
    if (DEBUG) {
        printf("Arbre apres la deuxieme passe:\n");
        printnoeud(racine, 0);
    }

    /*TROISIEME PASSE */
    /* A ce stade, on a un pseudocoarbre de racine racine,
       sans noeuds a un fils, et dont les etiquettes sont
       FEUIILLE ou UNKN. Il y a une pile des noeuds UNKN, stockes
       dans l'ordre de la postvisite d'un parcours en profondeur.
       On va s'en servir pour faire remonter les bords et les
       separateurs de bas en haut */

    taillepileinterne = indicepileinterne;
    for (indicepileinterne = 0; indicepileinterne < taillepileinterne;
         indicepileinterne++) {
        noeud *scanne;
        fils *nextfils;
        noeud *succourant;
        /* scanne est le noeud (pere) que l'on regarde;
           nextfils parcours sa liste de fils;
           courant est le fils actuellement examine et
           succourant=succ(courant) */
        noeud *debutnoeud;
        fils *debutfils;
        /*deux variables utilise pour la recherche de jumeaux:
           debut du bloc max */

        scanne = pileinterne[indicepileinterne];        /*he oui, la pile */
        nextfils = scanne->fils;        /*on commence au premier fils */
        do {
            /*la boucle chiante... cf mon memoire de DEA */
            courant = nextfils->pointe;
            /* bords */
            scanne->bg = min(scanne->bg, courant->bg);
            scanne->bd = max(scanne->bd, courant->bd);
            /*separateurs */
            scanne->ps = min(scanne->ps, courant->ps);
            if (scanne->fils->pointe != courant)
                /*ce n'est pas le premier fils */
                scanne->ps = min(scanne->ps, ps[(courant->bg) - 1]);
            scanne->ds = max(scanne->ds, courant->ds);
            if (scanne->lastfils->pointe != courant)
                /*ce n'est pas le dernier fils */
                scanne->ds = max(scanne->ds, ds[courant->bd]);

            nextfils = nextfils->suiv;
        }
        while (nextfils != NULL);


        if (DEBUG)
            printf("noeud %i-%i: ps=%i ds=%i", 1 + scanne->bg,
                   1 + scanne->bd, 1 + scanne->ps, 1 + scanne->ds);

        /* maintenant le test tout simple pour savoir si on a un module: */
        if (((scanne->ps) == (scanne->bg))
            && ((scanne->ds) == (scanne->bd))) {
            /*on a un module */
            scanne->type = MODULE;
            if (DEBUG)
                printf(" Module.\n");
        }
        else {
            scanne->type = ARTEFACT;
            if (DEBUG)
                printf(" artefact.\n");
        }
    }

    if (DEBUG) {
        printf("Arbre apres la troisieme passe:\n");
        printnoeud(racine, 0);
    }

    /* QUATRIEME PASSE */
    /* technique:on se contente de fusionner les artefacts a leurs parents
       ca se fait de bas en haut grace a pileinterne (toujours elle!) */

    for (indicepileinterne = 0; indicepileinterne < taillepileinterne;
         indicepileinterne++) {
        noeud *scanne;
        scanne = pileinterne[indicepileinterne];        /*he oui, la pile */
        if (scanne->type == ARTEFACT) {
            /*attention! La pile peut contenir des noeuds deconnectes */
            fusionne(scanne->pere, scanne);
            if (DEBUG)
                printf("Artefact elimine: %i-%i\n", 1 + scanne->bg,
                       1 + scanne->bd);
        }
    }
    if (DEBUG) {
        printf("Arbre apres la quatrieme passe:\n");
        printnoeud(racine, 0);
    }

    /* CINQIEME ET DERNIERE PASSE */
    /* on va typer les noeuds et extraire les fusions.
       comment on fait? Ben, la pile.... */
    for (indicepileinterne = 0; indicepileinterne < taillepileinterne;
         indicepileinterne++) {
        noeud *scanne;
        fils *nextfils;
        noeud *succourant;
        /* scanne est le noeud (pere) que l'on regarde;
           nextfils parcours sa liste de fils;
           courant est le fils actuellement examine et
           succourant=succ(courant) */
        noeud *debutnoeud;
        fils *debutfils;
        /*deux variables utilise pour la recherche de jumeaux:
           debut du bloc max */

        scanne = pileinterne[indicepileinterne];
        if (scanne->type != MODULE)
            continue;                /* je traite que les modules */

        nextfils = scanne->fils;        /*on commence au premier fils */
        while (1) {
            courant = nextfils->pointe;
            succourant = nextfils->suiv->pointe;
            if (ps[courant->bd] >= courant->bg
                && ds[courant->bd] <= succourant->bd) {
                /*Des jumeaux!! ->serie ou parallele! */
                /* on va determiner le bloc max de jumeaux consecutifs */
                debutnoeud = courant;
                debutfils = nextfils;
                while (ps[courant->bd] >= courant->bg &&
                       ds[courant->bd] <= succourant->bd &&
                       nextfils->suiv != NULL) {
                    nextfils = nextfils->suiv;
                    courant = nextfils->pointe;
                    if (nextfils->suiv == NULL)
                        break;
                    succourant = nextfils->suiv->pointe;
                }
                /*maintenant on connait la taille du bloc: il va de
                   debutnoeud a courant inclus,
                   et dans la liste des fils de scanne,
                   il s'etend de debutfils a nextfils inclus.
                   On extrait cette liste pour en faire les fils d'un
                   nouveau noeud... sauf si pas la peine! */
                if (debutfils == scanne->fils
                    && nextfils == scanne->lastfils) {
                    /* le noeud scanne etait serie ou parallele */
                  if ( adjii[debutnoeud->bd] !=0)
                        scanne->type = SERIE;
                    else
                        scanne->type = PARALLELE;
                }
                else {
                  if ( adjii[debutnoeud->bd]!=0)
                    /*seule cette ligne distingue G de G' !! */
                    {
                        nouveau = nouvnoeud(SERIE, scanne, -1, n);
                        if (DEBUG)
                            printf("Module serie extrait: %i-%i\n",
                                   1 + debutnoeud->bg, 1 + courant->bd);
                    }
                    else {
                        nouveau = nouvnoeud(PARALLELE, scanne, -1, n);
                        if (DEBUG)
                            printf("Module parallele extrait: %i-%i\n",
                                   1 + debutnoeud->bg, 1 + courant->bd);
                    }
                    extraire(scanne, nouveau, debutfils, nextfils);
                }
            }
            if (scanne->type == MODULE)
                scanne->type = PREMIER;
            if (nextfils->suiv == NULL || nextfils->suiv->suiv == NULL
                || nextfils->suiv->suiv->suiv == NULL)
                break;
            nextfils = nextfils->suiv;
        }
    }
    if (DEBUG) {
        printf("Arbre final:\n");
        printnoeud(racine, 0);
    }
    return racine;
}

void PrintG(graphe G)
/* affiche le graphe */
{
  int i,r; adj *a;
  r=0;
  for(i=0;i<G.n;i++)
    {
      printf("%i : ",i);
      a=G.G[i];
      while(a!=NULL)
        {
          printf("%i ", a->s);
          a=a->suiv;
        }
      printf("\n");
    }
}
void PrintGS(sommet **S, int n)
/* affiche le graphe trie selon S (celui utilise par algo2)*/
{
  int i; sadj *a;
  for(i=0;i<n;i++)
    {
      printf("%i : ",i);
      a=S[i]->adj;
      while(a!=NULL)
        {
          printf("%i ", a->pointe->place);
          a=a->suiv;
        }
      printf("\n");
    }
}

void PrintS2(sommet **S, int n)
     /* affiche la permutation factorisante */
{
  int i;
  printf(  "Place (nouvelle num) ");
  for(i=0;i<n;i++)
    printf("%3i ",S[i]->place);
  printf("\nNom (ancienne num) : ");
  for(i=0;i<n;i++)
    printf("%3i ",S[i]->nom);
  printf("\n");
}


/* la fonction principale; qui fait pas grand'chose....*/
noeud *decomposition_modulaire(graphe G)
{

    sommet **S;                        /* la permutation factorisante */
    noeud *Racine;                /* le futur arbre de decomposition */

    setbuf(stdout,NULL);

    S = algo1(G);               /* premiere partie: calcul
                                   de la permutation factorisante */

    TrierTous(S,G.n,Calculm(G));/* Trie les listes d'adjacence selon S
                                 */
    if(DEBUG)
      {
        PrintGS(S,G.n);
        PrintS2(S,G.n);
      }
    Racine = algo2(G, S);       /* deuxieme partie: calcul de l'arbre */
    return Racine;
}
