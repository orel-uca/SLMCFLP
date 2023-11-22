
#include <ilcplex/ilocplex.h> 
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <algorithm>
#include <vector>     
#include <iostream>
#include <fstream>
#define EPS 0.01

typedef IloArray<IloNumArray> FloatMatrix;
typedef IloArray<FloatMatrix> Float3Matrix;
typedef IloArray<IloExprArray> ExprMatrix;
typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<NumVarMatrix> NumVar3Matrix;

ILOSTLBEGIN

void SetParameters(IloCplex cplex)  
{
	cplex.setParam(IloCplex::TiLim, 200); 
	cplex.setOut(cplex.getEnv().getNullStream()); 
	
	cplex.setParam(IloCplex::FPHeur, 1);
	cplex.setParam(IloCplex::BndStrenInd, 1);
	cplex.setParam(IloCplex::MIPEmphasis, 4);
}



ILOMIPINFOCALLBACK7(IncumbentFeasibleCallback,
	IloCplex, cplex,
	IloBool, aborted,
	IloNum, tol,
	IloNum, lastIncumbent,
	IloNum, lastIncumbentTime,
	IloNum, timeStart,
	IloNum, timeLimit)  
{
	if (aborted == IloFalse && hasIncumbent() == IloFalse) { 
		IloNum timeUsed = cplex.getCplexTime() - timeStart;
		if (timeUsed > timeLimit + 20) { 
			getEnv().out() << "No feasible solution found in IncumbentFeasible (IncumbentFeasibleCallback)"
				<< timeUsed << endl;
			aborted = IloTrue;
			abort();
		}
	}
	if (!aborted) {
		if (hasIncumbent() &&
			fabs(lastIncumbent - getIncumbentObjValue()) > tol) {
			lastIncumbent = getIncumbentObjValue();
			lastIncumbentTime = cplex.getCplexTime();
		}
		else {
			if (hasIncumbent() && cplex.getCplexTime() - lastIncumbentTime > timeLimit) {
				getEnv().out() << "No improvement in the incumbent solution found (in IncumbentFeasible)"
					<< timeLimit << endl;
				aborted = IloTrue;
				abort();
			}
		}
	}
}

void main(int argc, char* argv[])
{
	char* datos = argv[1];
	char* output = argv[2];
	char* output2 = argv[3];
	
	

	IloEnv env;

	try
	{

		
		IloInt no, np, na, nc, DV_B=0;
		IloNum B, aux;
		IloNum tol = 0.1;
		IloNum tol_UBLB = 0.0001;
		IloNum tol_k0 = 0.35; 
		IloNum TL_kernel = 1800;
		float maxtime = 7200;
		IloNum checkcut = 0;
		IloNum UB = 0.0, LB = 0.0;
		IloNum tiempoFacKernel = 10.0;

		ifstream entrada(datos);
		if (!datos) {
			cout << "Cannot open file.\n";
			system("pause");
			throw(-1);
		}
		entrada.ignore(3);
		entrada >> no; 
		entrada.ignore(5, ':');
		entrada >> np; 
		entrada.ignore(8, ':');
		entrada >> na; 
		entrada.ignore(11, ':');
		entrada >> nc; 
		entrada.ignore(14, ':');
		entrada >> B; 
		env.out() << datos << endl;
		env.out() << no << " " << np << " " << na << " " << nc << " " << B << endl;

		IloNumArray Rfa(env, np); 
		IloNumArray rac(env, np); 
		IloNumArray m(env, no); 
		IloNumArray bj(env, na);
		FloatMatrix c(env, na);
		for (int j = 0; j < na; j++) c[j] = IloNumArray(env, np);
		FloatMatrix w(env, nc);
		for (int i = 0; i < nc; i++) w[i] = IloNumArray(env, no);
		FloatMatrix d(env, no + na + nc);
		for (int i = 0; i < no + na + nc; i++) d[i] = IloNumArray(env, no + na + nc);

		entrada.ignore(22, '[');
		for (int i = 0; i < np; i++) {
			entrada >> rac[i];
			entrada >> Rfa[i];
		}
		entrada.ignore(15, '[');
		for (int s = 0; s < no; s++) {
			entrada >> aux;
			entrada >> aux; 
			entrada >> m[s]; 
			m[s] -= 1; 
		}

		entrada.ignore(15, '[');
		for (int j = 0; j < na; j++) {
			entrada >> aux; 
			entrada >> aux; 
			entrada >> bj[j]; 
			for (int t = 0; t < np; t++) {
				entrada >> c[j][t];
			}
		}
		entrada.ignore(15, '[');
		for (int i = 0; i < nc; i++) {
			entrada >> aux;
			entrada >> aux; 
			for (int s = 0; s < no; s++) {
				entrada >> w[i][s];
			}
		}
		entrada.ignore(15, '[');
		for (int i = 0; i < no + na + nc; i++) {
			for (int j = 0; j < no + na + nc; j++) entrada >> d[i][j];
		}
		entrada.close();
		env.out() << "Los archivos fueron cargados correctamente" << endl;


		ofstream salida(output, ios::app);
		ofstream salida2(output2, ios::app);



		//----------------------------//
		IloEnv env1;
		IloNum num_aux = 0;
		FloatMatrix exists_y(env1, nc); 
		FloatMatrix exists_y_K(env1, nc);
		for (int i = 0; i < nc; i++) {
			exists_y[i] = IloNumArray(env1, no);
			exists_y_K[i] = IloNumArray(env1, no);
			for (int s = 0; s < no; s++) {
				exists_y[i][s] = 0;
				aux = 0;
				for (int j = 0; j < na; j++) {
					if ((d[s][no + j] <= Rfa[m[s]]) && (d[no + j][no + na + i] <= rac[m[s]])) aux = 1;
				}
				if (aux == 1) {
					exists_y[i][s] = 1; 
					num_aux +=1;
				}
			}
		}

		IloNum num_y = num_aux;
		num_aux = 0;

		FloatMatrix exists_x(env1, na); 
		FloatMatrix exists_x_K(env1, na);
		for (int j = 0; j < na; j++) {
			exists_x[j] = IloNumArray(env1, np);
			exists_x_K[j] = IloNumArray(env1, np);
			for (int mp = 0; mp < np; mp++) {
				exists_x[j][mp] = 0;
				aux = 0;
				for (int s = 0; s < no; s++) {
					if ((d[(no + j)][s] <= Rfa[mp]) && (mp == m[s])) aux += 1;
				}
				if (aux >= 0.5) {
					exists_x[j][mp] = 1; 
					num_aux += 1;
				}
			}
		}

		IloNum num_x = num_aux;

		env1.out()<<num_y<<" " << num_x<<endl;
	

		////------------------------------------------------------------------------- //  
		time_t start, end, startk, endk, startm, endm;
		double total_time, time_kernel, time_m=0.0;
		time_t now = time(0);
		char* dt = ctime(&now);
		start = clock();
		//--------------------------------------------------------------------------------//
		
		///----------------------------------------------------------------------------------///
		
		IloModel modelo_k0(env1);

		NumVarMatrix y_k0(env1, nc); 
		for (int i = 0; i < nc; i++) {
			y_k0[i] = IloNumVarArray(env1, no);
			for (int s = 0; s < no; s++) {
				
				y_k0[i][s] = IloNumVar(env1, 0.0, 1.0, ILOFLOAT); 
			}
		}

		NumVarMatrix x_k0(env1, na); 
		for (int j = 0; j < na; j++) {
			x_k0[j] = IloNumVarArray(env1, np);
			for (int mp = 0; mp < np; mp++) {
				
				x_k0[j][mp] = IloNumVar(env1, 0.0, 1.0, ILOFLOAT); 
			}
		}

		NumVarMatrix u_k0(env1, na); 
		for (int j = 0; j < na; j++) {
			u_k0[j] = IloNumVarArray(env1, bj[j]);
			for (int mp = 0; mp < bj[j]; mp++) {
				
				u_k0[j][mp] = IloNumVar(env1, 0.0, 1.0, ILOFLOAT); 
			}
		}




		IloExpr fo_k0(env1);

		for (int i = 0; i < nc; i++) {
			for (int s = 0; s < no; s++) {
				if (exists_y[i][s] == 1) fo_k0 += w[i][s] * y_k0[i][s];
			}
		}
		IloObjective obj_k0 = IloMaximize(env1, fo_k0);
		modelo_k0.add(obj_k0);

		
		for (int i = 0; i < nc; i++) {
			for (int s = 0; s < no; s++) {
				if (exists_y[i][s] == 1) {
					IloExpr convsum(env1);
					for (int j = 0; j < na; j++) {
						if ((exists_x[j][m[s]] == 1) && (d[s][no + j] <= Rfa[m[s]]) && (d[no + j][no + na + i] <= rac[m[s]])) convsum += x_k0[j][m[s]];
					}
					modelo_k0.add(y_k0[i][s] <= convsum);
					convsum.end();
				}
			}
		}

		
		for (int i = 0; i < nc; i++) {
			for (int mp = 0; mp < np; mp++) {
				IloExpr con_maxprod(env1);
				for (int s = 0; s < no; s++) {
					if ((m[s] == mp) & (exists_y[i][s] == 1)) con_maxprod += y_k0[i][s];
				}
				modelo_k0.add(con_maxprod <= 1);
				con_maxprod.end();
			}
		}

		
		IloExpr con_B_k0(env1);
		
		for (int j = 0; j < na; j++) {
			
			IloExpr con_numprod(env1);
			
			IloExpr con_tam(env1);
			for (int mp = 0; mp < np; mp++) {
				if (exists_x[j][mp] == 1) con_numprod += x_k0[j][mp];
			}
			for (int mp = 0; mp < bj[j]; mp++) {
				con_numprod -= (mp + 1) * u_k0[j][mp];
				con_tam += u_k0[j][mp];
				con_B_k0 += c[j][mp] * u_k0[j][mp];
			}
			modelo_k0.add(con_numprod == 0);
			modelo_k0.add(con_tam <= 1);
			con_numprod.end();
			con_tam.end();
		}
		modelo_k0.add(con_B_k0 <= B);
		con_B_k0.end();

		IloNumArray auxJ(env1, na);
		for (int j = 0; j < na; j++) auxJ[j] = 0;

		for (int i = 0; i < nc; i++) {
			for (int t = 0; t < np; t++) {
				for (int s = 0; s < no; s++) {
					if ((exists_y[i][s] == 1) & (m[s] == t)) {
						IloExpr new2(env1);

						for (int ss = 0; ss < no; ss++) {
							if ((m[ss] == t) & (exists_y[i][ss] == 1) & (w[i][ss] >= w[i][s])) {
								new2 += y_k0[i][ss];
								for (int j = 0; j < na; j++) {
									if ((d[ss][no + j] <= Rfa[1]) & (d[no + j][no + na + i] <= rac[1]) & (exists_x[j][t] == 1)) {
									
										auxJ[j] = 1;
									}
								}
							}
						}

					
						for (int j = 0; j < na; j++) {
							if (auxJ[j] == 1) new2 -= x_k0[j][t];
						}
						for (int j = 0; j < na; j++) auxJ[j] = 0;

						if (new2.isConstant() == 0) {
							
							checkcut = 1;
							modelo_k0.add(new2 <= 0);
							new2.end();
						}
					}
				}
			}
		}

		
		IloCplex cplex_k0(modelo_k0);

		cplex_k0.setOut(env.getNullStream());
		cplex_k0.setParam(IloCplex::TiLim, maxtime);
		

		startk = clock();
		cplex_k0.solve();
		endk = clock();
		time_kernel = (double)(endk - startk) / (double)CLOCKS_PER_SEC;

		
		

		salida << datos<< " "<< cplex_k0.getStatus() << " ";
		salida2 << datos << " " << cplex_k0.getStatus() << " ";
		
		if (cplex_k0.getStatus() == IloAlgorithm::Status::Feasible || cplex_k0.getStatus() == IloAlgorithm::Status::Optimal) {
			salida << cplex_k0.getObjValue() << " ";
			salida2<< cplex_k0.getObjValue() << " "<<endl;
		
		}
		else {
			cout << "Unfeasible LP.\n";
			system("pause");
			throw(-1);
		}

		UB = cplex_k0.getObjValue();

		//----------------------
		
		//----------------------
		vector<int> K0_al;
		vector<float> rcost_al(na,-10000);
		
		for (int j = 0; j < na; j++) {
			for (int mp = 0; mp < np; mp++) {
				if (exists_x[j][mp] == 1) {
					if (cplex_k0.getValue(x_k0[j][mp] >= tol_k0)) {
						
						K0_al.push_back(j);
						rcost_al[j] = mp * tol_k0 + 1; 
						break; 
					}
					else {
						if (cplex_k0.getValue(x_k0[j][mp] >= EPS)) rcost_al[j] = IloMax(cplex_k0.getValue(x_k0[j][mp]), rcost_al[j]);
						else rcost_al[j]=IloMax(cplex_k0.getReducedCost(x_k0[j][mp]), rcost_al[j]);
					}
				}
			}
		}

		cplex_k0.end();
		modelo_k0.end();
		obj_k0.end();
		fo_k0.end();

		
		
		
		

		vector<size_t> idx_2(rcost_al.size());
		for (size_t i = 0; i != idx_2.size(); ++i) idx_2[i] = i;
		sort(idx_2.begin(), idx_2.end(), [&rcost_al](size_t i1, size_t i2) {return rcost_al[i1] > rcost_al[i2]; });
		rcost_al.end();
		
		idx_2.erase(idx_2.begin(), idx_2.begin()+K0_al.size()); 

		
		
		vector<int> K0;
		vector<int>  v_aux, kernel;
		for (int j = 0; j < K0_al.size(); j++) {
			for (int mp = 0; mp < np; mp++) {
				if (exists_x[K0_al[j]][mp] == 1) {
					K0.push_back(K0_al[j]*na+mp);
				}
			}
		}
		kernel.insert(kernel.begin(), K0.begin(), K0.end());
		IloNum tam_kernel;  
		IloNum minnumki;
		if ((na - idx_2.size()) > 0.3 * na) {
			tam_kernel = 2;
			minnumki = 2;
		}
		else {
			tam_kernel = 4;
			minnumki = IloMax(ceil((0.1 * na - (na - idx_2.size())) / 4),2);  
		}
		IloNum tam_idx = idx_2.size();
		
		
		
		
		vector<int> aux2(kernel.size(), 0); 
		
		IloNum numki = 0;
		int tamaux;
		IloNum ant_factible = 0;
		IloNum tiempo_ant = tiempoFacKernel+1; 
		IloNum aux_sum = 0;
		IloNum nuevas;
		IloNum num_infactible = 0; 

		salida2 << " min iter "<< minnumki << " tam_kernel "<< tam_kernel << endl; 
		//------------ 
		//-----------
				
		NumVarMatrix y_k(env1, nc);
		for (int i = 0; i < nc; i++) {
			y_k[i] = IloNumVarArray(env1, no);
			for (int s = 0; s < no; s++) {
				if (exists_y[i][s] == 1) {
					y_k[i][s] = IloNumVar(env1, 0.0, 1.0, ILOBOOL); 
				}
			}
		}

		NumVarMatrix x_k(env1, na);
		NumVarMatrix u_k(env1, na);  
		for (int j = 0; j < na; j++) {
			aux_sum = IloSum(exists_x[j]);
			if (aux_sum >= 1) {
				x_k[j] = IloNumVarArray(env1, np);
				u_k[j] = IloNumVarArray(env1, bj[j]);
				for (int mp = 0; mp < np; mp++) {
					
					if (exists_x[j][mp] == 1) {
						
						x_k[j][mp] = IloNumVar(env1, 0.0, 1.0, ILOBOOL); 
						
					}
					if (mp < IloMin(bj[j], aux_sum)) u_k[j][mp] = IloNumVar(env1, 0.0, 1.0, ILOBOOL); 
				}
			}
		}

		//----------
		
		//--------------
		do {
			if (((double)(clock() - start)/ (double)CLOCKS_PER_SEC) >= TL_kernel) break;
			numki += 1;
			
			if (abs(UB - LB)/(0.0000000001+ abs(UB)) < tol_UBLB) { 
				
				salida2 << " UB==LB" << endl;
					break;
			}
			if (tiempo_ant < tiempoFacKernel) { tam_kernel +=2 ; } 
			
			
			if (numki >= 2) {
				

				aux_sum = idx_2.size();
				if (aux_sum >= 1) {
					nuevas = IloMin(tam_kernel, aux_sum);
					for (int j = 0; j < nuevas; j++) {
						for (int mp = 0; mp < np; mp++) {
							if (exists_x[idx_2[j]][mp] == 1) {
								v_aux.push_back(idx_2[j] * na + mp);
							}
						}
					}
					idx_2.erase(idx_2.begin(), idx_2.begin() + nuevas);
					kernel.clear();
					kernel.insert(kernel.begin(), K0.begin(), K0.end());
					kernel.insert(kernel.end(), v_aux.begin(), v_aux.end());
				}
				else {
					salida2 << "No hay más candidatos que probar" << endl;
					
					break;
				}
				
				

			}
			////-----------------

			IloModel modeloK(env1);
			salida2 << "Iter " << numki << " ";
			//---------------
			//---------------


			
			for (int j = 0; j < na; j++) {
				for (int mp = 0; mp < np; mp++) {
					exists_x_K[j][mp] = 0;
				}
			}
			for (int j = 0; j < kernel.size(); j++) exists_x_K[(kernel[j] / na)][(kernel[j] % na)] = 1;

			
			for (int i = 0; i < nc; i++) {
				for (int s = 0; s < no; s++) {
					exists_y_K[i][s] = 0;
					aux = 0;
					for (int j = 0; j < kernel.size(); j++) { 
						if (kernel[j] % na==m[s]) { 
							if ((d[s][no + (kernel[j] / na)] <= Rfa[m[s]]) && (d[no + (kernel[j] / na)][no + na + i] <= rac[m[s]])) {
								aux = 1;
							}
						}
					}
					if (aux == 1) {
						exists_y_K[i][s] = 1; 
						num_aux += 1;
					}
				}
			}

			IloExpr fo_k(env1);

			for (int i = 0; i < nc; i++) {
				for (int s = 0; s < no; s++) {
					if (exists_y_K[i][s] == 1) fo_k += w[i][s] * y_k[i][s];
				}
			}
			IloObjective obj_k = IloMaximize(env1, fo_k);
			modeloK.add(obj_k);

			
			for (int i = 0; i < nc; i++) {
				for (int s = 0; s < no; s++) {
					if (exists_y_K[i][s] == 1) {
						IloExpr convsum(env1);
						for (int j = 0; j < na; j++) {
							if ((exists_x_K[j][m[s]] == 1) && (d[s][no + j] <= Rfa[m[s]]) && (d[no + j][no + na + i] <= rac[m[s]])) convsum += x_k[j][m[s]];
						}
						modeloK.add(y_k[i][s] <= convsum);
						convsum.end();
					}
				}
			}

			
			
			for (int i = 0; i < nc; i++) {
				if (IloSum(exists_y_K[i]) >= 1) {
					for (int mp = 0; mp < np; mp++) {
						aux_sum = 0;
						IloExpr con_maxprod(env1);
						for (int s = 0; s < no; s++) {
							if ((m[s] == mp) & (exists_y_K[i][s] == 1)) {
								con_maxprod += y_k[i][s];
								aux_sum += 1;
							}
						}
						if(aux_sum>=2) modeloK.add(con_maxprod <= 1);
						con_maxprod.end();
					}
				}
			}

			
			IloExpr con_B_K(env1);
			
			for (int j = 0; j < na; j++) {
				aux_sum = IloSum(exists_x_K[j]); 
				if (aux_sum>=1) {
					
					IloExpr con_numprod(env1);
					
					IloExpr con_tam(env1);
					for (int mp = 0; mp < np; mp++) {
						if (exists_x_K[j][mp] == 1) con_numprod += x_k[j][mp];
					}
					for (int mp = 0; mp < IloMin(bj[j], aux_sum ); mp++) {
						con_numprod -= (mp + 1) * u_k[j][mp];
						con_tam += u_k[j][mp];
						con_B_K += c[j][mp] * u_k[j][mp];
					}
					modeloK.add(con_numprod == 0);
					if(aux_sum >= 2) modeloK.add(con_tam <= 1);
					con_numprod.end();
					con_tam.end();
				}
			}
			modeloK.add(con_B_K <= B);
			con_B_K.end();

			
			if (numki >= 2) {
				
				IloExpr constAtLeastOne(env1);
				for (int j = 0; j < v_aux.size(); j++) {
					constAtLeastOne += x_k[v_aux[j]/na][v_aux[j] % na];
				}
				if (ant_factible == 1) {
					for (int j = 0; j < K0.size(); j++) {
						if (aux2[j] == 1) { 
							constAtLeastOne += x_k[K0[j] / na][K0[j] % na];
						}
					}
				}
				modeloK.add(constAtLeastOne >= 1);
				constAtLeastOne.end();
			}

			

			IloNumArray auxJ(env1, na);
			for (int j = 0; j < na; j++) auxJ[j] = 0;

			for (int i = 0; i < nc; i++) {
				for (int t = 0; t < np; t++) {
					for (int s = 0; s < no; s++) {
						if ((exists_y_K[i][s] == 1) & (m[s] == t)) {
							IloExpr new2(env1);

							for (int ss = 0; ss < no; ss++) {
								if ((m[ss] == t) & (exists_y_K[i][ss] == 1) & (w[i][ss] >= w[i][s])) {
									new2 += y_k[i][ss];
									for (int j = 0; j < na; j++) {
										if ((d[ss][no + j] <= Rfa[1]) & (d[no + j][no + na + i] <= rac[1]) & (exists_x_K[j][t] == 1)) {
											
											auxJ[j] = 1;
										}
									}
								}
							}

							
							for (int j = 0; j < na; j++) {
								if (auxJ[j] == 1) new2 -= x_k[j][t];
							}
							for (int j = 0; j < na; j++) auxJ[j] = 0;

							if (new2.isConstant() == 0) {
								
								checkcut = 1;
								modeloK.add(new2 <= 0);
								new2.end();
							}
						}
					}
				}
			}

			IloCplex cplexK(modeloK);
			SetParameters(cplexK);
			cplexK.setParam(IloCplex::Param::MIP::Tolerances::LowerCutoff, LB);
			cplexK.use(IncumbentFeasibleCallback(env1, cplexK, IloFalse, 0.001, IloInfinity, cplexK.getCplexTime(), cplexK.getCplexTime(), 100)); 
		

			startk = clock();
			cplexK.solve();
			endk = clock();
			time_kernel += (double)(endk - startk) / (double)CLOCKS_PER_SEC;
			tiempo_ant = (double)(endk - startk) / (double)CLOCKS_PER_SEC;

			

			if (cplexK.getStatus() == IloAlgorithm::Infeasible || cplexK.getStatus() == IloAlgorithm::Unknown) {
				num_infactible += 1;
				salida2 << cplexK.getStatus() << " num_bin_var " << cplexK.getNbinVars() << " ";
				cplexK.end();
				modeloK.end();
				obj_k.end();
				fo_k.end();
				if (num_infactible == 3) break;
			}
			else {
				if (cplexK.getStatus() == IloAlgorithm::Feasible) { ant_factible = 1; }
				num_infactible = 0;
				LB = cplexK.getObjValue();
				salida2 << cplexK.getStatus() <<" Solución " << LB << " num_bin_var " << cplexK.getNbinVars() << " ";
				

				if (numki >= 2) {
					for (int j = 0; j < K0.size(); j++) {
						if (cplexK.getValue(x_k[K0[j] / na][K0[j] % na]) < tol) {
							aux2[j] += 1;
							
							if (aux2[j] == 2) { 
								K0.erase(K0.begin() + j);
								aux2.erase(aux2.begin() + j);
								j = j - 1; 
							}
						}
						else aux2[j] = 0; 
					}
					for (int j = 0; j < v_aux.size(); j++) { 
						if (cplexK.getValue(x_k[v_aux[j] / na][v_aux[j] % na]) > tol) {
							K0.push_back(v_aux[j]);
							aux2.push_back(0);
						}
					}
					v_aux.clear();
				}

				
				startm = clock();
				cplexK.end();
				modeloK.end();
				obj_k.end();
				fo_k.end();
				endm = clock();
				time_m += (double)(endm - startm) / (double)CLOCKS_PER_SEC;
			}

			
			salida2 << "tiempo " << ((double)(clock() - start) / (double)CLOCKS_PER_SEC) << endl; 
			
		}while ((numki < minnumki) || (((double)(clock() - start) / (double)CLOCKS_PER_SEC ) <100) ); 
		

		end = clock();
		
		total_time = (double)(end - start) / (double)CLOCKS_PER_SEC;
		salida << LB<<" "<< total_time <<" " << time_kernel<< " "<< numki<< " ";
		salida2 << "LB " <<LB << " Tiempo total " << total_time << " Tiempo resolución subproblemas " << time_kernel << " Tiempo liberar memoria " << time_m<<endl;
		
        salida << endl;
		env1.end();
		
	

	}


	catch (IloException& ex)
	{
		cerr << "Error Cplex: " << ex << endl;
	}

	catch (...)
	{
		cerr << "Error Cpp" << endl;
	}

	env.end();

	
}
