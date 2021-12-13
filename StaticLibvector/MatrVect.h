#pragma once
#include<iostream>
#include<assert.h>
#include<math.h>
#include "pch.h"
using namespace std;
namespace mat {
	template<class T>
	class Vector;
	template<class T>
	class Matrice {
	private:
		unsigned int nb_col,nb_ligne;
		T** matr;
	public:
		Matrice(unsigned int nbl, unsigned int nbc)
		{
			this->nb_col = nbc;
			this->nb_ligne = nbl;
			this->matr = new int* [nbl];
			for (int i = 0; i < nbl; i++) this->matr[i] = new int[nbc];
			//assert(nbc > 0 && nbl > 0);
			try
			{
				if (nbc <= 0 || nbl <= 0) throw " Matrice dimension erre";
				for (int i = 0; i < nbl; i++)
				{
					for (int j = 0; j < nbc; j++) this->matr[i][j] = 0;
				}
			}
			catch (char* e)
			{
				cout << e << endl;
			}

		}
		Matrice(Matrice& m) = delete;
		~Matrice()
		{
			if (this->matr)
			{
				for (int i = 0; i < this->nb_ligne; i++)
				{
					delete[] this->matr[i];
					this->matr[i] = 0;
				}
				delete[] this->matr;
				this->matr = NULL;
			}
		}

		void remplissage(T val)
		{
			for (int i = 0; i < this->nb_ligne; i++)
			{
				for (int j = 0; j < this->nb_col; j++) this->matr[i][j] = val;
			}
		}
		Matrice& operator=(const Matrice& m)= delete;
		Matrice* operator+(const Matrice& m)const
		{
			//assert(this->nb_col == m.nb_col && this->nb_ligne == m.nb_ligne);
			try
			{
				if (this->nb_col != m.nb_col || this->nb_ligne != m.nb_ligne) throw " impossible d'aditionner 2 matrices qui n'ont pas la meme dimension";
				Matrice* m_addit = new Matrice(this->nb_ligne, this->nb_col);
				for (int i = 0; i < this->nb_ligne; i++)
					for (int j = 0; j < this->nb_ligne; j++)m_addit->matr[i][j] = this->matr[i][j] + m.matr[i][j];
				return m_addit;
			}
			catch (char* e)
			{
				cout << e << endl;
			}

			// TODO: insérer une instruction return ici
		}
		Matrice* operator-(const Matrice& m)const
		{
			try
			{
				if (this->nb_col != m.nb_col || this->nb_ligne != m.nb_ligne) throw " impossible de soustracter 2 matrices qui n'ont pas la meme dimension";
				Matrice* m_sous = new Matrice(this->nb_ligne, this->nb_col);
				for (int i = 0; i < this->nb_ligne; i++)
					for (int j = 0; j < this->nb_ligne; j++)m_sous->matr[i][j] = this->matr[i][j] - m.matr[i][j];
				return m_sous;
			}
			catch (char* e)
			{
				cout << e << endl;
			}

			// TODO: insérer une instruction return ici
		}
		Matrice* operator*(const Matrice& m)const
		{
			try
			{
				if (this->nb_col != m.nb_ligne) throw " impossible de multiplier ces 2 matrices dimension error";
				Matrice* m_mult = new Matrice(this->nb_ligne, m.nb_col);
				for (int i = 0; i < this->nb_ligne; i++)
					for (int j = 0; j < m.nb_col; j++)
						for (int k = 0; k < m.nb_ligne; k++)m_mult->matr[i][j] += this->matr[i][k] * m.matr[k][j];
				return m_mult;
			}
			catch (char* e)
			{
				cout << e << endl;
			}
			// TODO: insérer une instruction return ici
		}
		int* operator[](unsigned int ind1)
		{

			return this->matr[ind1];

		}
		//Matrice* operator[](unsigned int ind2);
		Vector<T>* matr_to_vect()
		{
			assert(this->nb_ligne == 1 && this->nb_col == 3);
			Vector<T>* v = new Vector<T>(this->matr[0][0], this->matr[0][1], this->matr[0][2]);
			return v;
		}
		void print() const
		{
			for (int i = 0; i < this->nb_ligne; i++)
			{
				for (int j = 0; j < this->nb_col; j++) cout << "matrice[" << i << "]" << "[" << j << "]=" << this->matr[i][j] << "\t";
				cout << endl;
			}
		}

	};
	template<class T>
	class Vector
	{
	private:
		T coor[3];
	public:
		Vector(T x = 0.0, T y = 0.0, T z = 0.0)
		{
			this->coor[0] = x;
			this->coor[1] = y;
			this->coor[2] = z;
		}
		Vector(const Vector&) = delete;
		Vector& operator=(const Vector&) = delete;
		void print() const
		{
			for (int i = 0; i < 3; i++)
			{
				cout << "coor[" << i << "]=" << this->at(i) << endl;
			}
			//cout << "x=" << this->at(0)  << "\ty=" << this->at(1) << "\tz=" << this->at(2);
		}
		T at(unsigned int indice) const
		{
			assert(indice >= 0 && indice < 3);//debug
			try //release
			{
				if (indice < 0 || indice >= 3) throw "rang erre";
				return this->coor[indice];
			}
			catch (char* e)
			{
				cout << e << endl;
			}
		}

		T& operator[](unsigned int indice)
		{
			assert(indice >= 0 && indice < 3);//debug
			try //release
			{
				if (indice < 0 || indice >= 3) throw "rang erre";
				return this->coor[indice];
			}
			catch (char* e)
			{
				cout << e << endl;
			}
		}
		bool operator==(const Vector& v)const
		{
			bool res = true;
			for (int i = 0; i < 3; i++)
			{
				res = res && (this->at(i) == v.at(i));
			}
			return res;
		}
		bool operator!=(const Vector& v)const
		{
			bool res = !(*this == v);
			return res;
		}

		T module()const
		{
			return(sqrt(pow(this->at(0), 2) + pow(this->at(1), 2) + pow(this->at(2), 2)));
		}
		T operator*(const Vector& v)const  //le produit scalaire
		{
			return(this->at(0) * v.at(0) + this->at(1) * v.at(1) + this->at(2) * v.at(2));
		}
		void operator*(float val)			//le produit d'un vecteur et une valeur
		{
			this->coor[0] *= val;
			this->coor[1] *= val;
			this->coor[2] *= val;
		}
		Vector* operator^(const Vector& v)const//le produit vectoriel
		{
			Vector* prod = new Vector();
			prod->coor[0] = this->at(1) * v.at(2) - this->at(2) * v.at(1);
			prod->coor[1] = this->at(0) * v.at(2) - this->at(2) * v.at(0);
			prod->coor[2] = this->at(0) * v.at(1) - this->at(1) * v.at(0);
			return prod;
		}
		Vector* operator+(const Vector& v)const
		{
			Vector* add = new Vector();
			add->coor[0] = this->at(0) + v.at(0);
			add->coor[1] = this->at(1) + v.at(1);
			add->coor[2] = this->at(2) + v.at(2);
			return add;
		}
		Vector* operator-(const Vector& v)const
		{
			Vector* sous = new Vector();
			sous->coor[0] = this->at(0) - v.at(0);
			sous->coor[1] = this->at(1) - v.at(1);
			sous->coor[2] = this->at(2) - v.at(2);
			return sous;
		}
		void operator/(T val)
		{
			this->coor[0] /= val;
			this->coor[1] /= val;
			this->coor[2] /= val;
		}
		Matrice<T>* vect_to_matr()
		{
			Matrice<T>* m = new Matrice<T>(1, 3);
			(*m)[0][0] = this->coor[0];
			(*m)[0][1] = this->coor[1];
			(*m)[0][2] = this->coor[2];
			return m;
		}
		
	};
};
