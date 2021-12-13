// uservector.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
//

#include <iostream>
//#include"vector.h"
#include"MatrVect.h"
using namespace std;
using namespace mat;
int main()
{
    Matrice<int> A(2, 2);
    //Matrice* B=(1, 3);
    Matrice<int> C(1, 3);
    A.remplissage(5);
    A.print();
    A[0][1] = 9;
    A.print();
    Vector<int> v1(1, 2, 3);
    v1.print();
    Matrice<int>* B = v1.vect_to_matr();
    B->print();
    C.remplissage(6);
    C.print();
    Vector<int>* v2 = C.matr_to_vect();
    v2->print();
    /*B.remplissage(3);
    B.print();
    C.remplissage(2);
    C.print();
    Matrice* add =new Matrice(2,2);
    add = A + B;
    cout << "addition:\t" << endl;
    add->print();
   Matrice* sous=new Matrice(2, 2);
    sous=A - B;
    cout << "soustraction:\t"<< endl;
    sous->print();
    Matrice* mul=new Matrice(2, 3);
    mul=B * C;
    cout << "multiplication:"<< endl;
    mul->print();*/

    /*Vector v1(1, 2, 3);
    v1.print();
    v1[1] = 8;
    v1.print();
    Vector a(2, 5, 6);
    Vector b(1, 2, 3);
    a.print();
    b.print();
    Vector* tmp = new Vector();
    tmp = a + b;
    tmp->print();
    tmp = a - b;
    tmp->print();
    tmp = a ^ b;
    tmp->print();
    float module = a.module();
    cout << "\n le module=" << module << endl;
    a * 2;
    a.print();
    a / 2;
    a.print();
    delete tmp;*/
}

// Exécuter le programme : Ctrl+F5 ou menu Déboguer > Exécuter sans débogage
// Déboguer le programme : F5 ou menu Déboguer > Démarrer le débogage

// Astuces pour bien démarrer : 
//   1. Utilisez la fenêtre Explorateur de solutions pour ajouter des fichiers et les gérer.
//   2. Utilisez la fenêtre Team Explorer pour vous connecter au contrôle de code source.
//   3. Utilisez la fenêtre Sortie pour voir la sortie de la génération et d'autres messages.
//   4. Utilisez la fenêtre Liste d'erreurs pour voir les erreurs.
//   5. Accédez à Projet > Ajouter un nouvel élément pour créer des fichiers de code, ou à Projet > Ajouter un élément existant pour ajouter des fichiers de code existants au projet.
//   6. Pour rouvrir ce projet plus tard, accédez à Fichier > Ouvrir > Projet et sélectionnez le fichier .sln.
