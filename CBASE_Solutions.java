/************************************************************           
 * MutPanning 												*
 * 															*   
 * Author:		Felix Dietlein								*   
 *															*   
 * Copyright:	(C) 2019 									*   
 *															*   
 * License:		BSD-3-Clause open source license			*   
 *															*   
 * Summary: This script optimizes the parameters of the 	*
 * CBASE model to model the background distribution of		*
 * synonymous mutations across the exome. For this purpose,	*
 * we use a quasi-Newton approach (Limited-Memory BFGS) for	*
 * the optimization. For each model we start with 150 random*   
 * initial choices for the parameters, perform L-BFGS and 	*
 * select the best solution.								*
 *************************************************************/

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import org.apache.commons.math3.special.Gamma;
import jdistlib.math.Bessel;

public class CBASE_Solutions {
	static String[] index_header_samples={"ID","Sample","Cohort"};
	
	/*
	 * argument0 root file
	 * arugment1 sample file
	 */
	
	static double[] log_precompute=new double[1000000];
	static double[] gamma_precompute=new double[1000000];
	static double[] gamma_delta_precompute=new double[1000000];
	static double[] exp_precompute=new double[2000000];
//	static double[][] kv_precompute=new double[1000][4000];
//	static double[][] kv_precompute_diff1=new double[1000][4000];
//	static double[][] kv_precompute_diff2=new double[1000][4000];
	
	
	public static void main (String[] args){
		/*
		for (double i=0.1;i<=10;i+=0.1){
			System.out.print(i);
			for (double j=0.1;j<=10;j+=0.1){
				System.out.print("	"+besselk_diff2_short(i,j));
			}
			System.out.println();
		}
		System.exit(0);
		*/
		
		for (int i=0;i<gamma_precompute.length;i++){
			log_precompute[i]=Math.log((double)(i)/1000.0);
			gamma_precompute[i]=Gamma.logGamma((double)(i)/1000.0);
			gamma_delta_precompute[i]=(Gamma.logGamma(0.0001+(double)(i)/1000.0)-Gamma.logGamma(-0.0001+(double)(i)/1000.0))/(2*0.0001);
		}
		
		for (int i=0;i<exp_precompute.length;i++){
			exp_precompute[i]=Math.exp((double)(i)/1000.0-1000);
		}
	
		/*
		for (int i=0;i<kv_precompute.length;i++){
			for (int j=0;j<kv_precompute[i].length;j++){
				kv_precompute[i][j]=Bessel.k((double)(i)/100.0,(double)(j)/100.0-20,false);
				kv_precompute_diff1[i][j]=(Bessel.k((double)(i)/100.0,(double)(j)/100.0-20+0.0001,false)-Bessel.k((double)(i)/100.0,(double)(j)/100.0-20-0.0001,false))/(2*0.0001);
				kv_precompute_diff2[i][j]=(Bessel.k((double)(i)/100.0+0.0001,(double)(j)/100.0-20,false)-Bessel.k((double)(i)/100.0-0.0001,(double)(j)/100.0-20,false))/(2*0.0001);
				
			}
		}*/
		
		//Determine all entity names, as clustering is performed separately for each Step this is needed to coordinate the order
		String[] entity_name=new String[0];
		try{
			FileInputStream in=new FileInputStream(args[1]);
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header=index_header(input.readLine().split("	"),index_header_samples);
			String s="";
			ArrayList<String> aa =new ArrayList<String>();
			while((s=input.readLine())!=null){
				String e=s.split("	")[index_header[2]];
				if(!contains(e,aa)){
					aa.add(e);
				}
			}
			input.close();
			Collections.sort(aa);
			//aa.add("PanCancer");
			entity_name=new String[aa.size()];
			for (int i=0;i<aa.size();i++){
				entity_name[i]=aa.get(i);
			}
		}
		catch(Exception e){
			System.out.println(e);
		}
		
		Comparator<int[]> comp=(int[] a, int[] b)->{
			return new Integer(a[0]).compareTo(b[0]);
		};
		
		if(!new File(args[0]+"CBASE/Parameters_Summary/").exists()){
			new File(args[0]+"CBASE/Parameters_Summary/").mkdir();
		}
		
		for (int k=0;k<entity_name.length;k++){
			//if(k<=8){
			//	continue;
			//}
			System.out.println(entity_name[k]);
			ArrayList<Integer> counts=new ArrayList<Integer>();
			try{
				FileInputStream in=new FileInputStream(args[0]+"CBASE/Counts/CountSilent"+entity_name[k]+".txt");
				DataInputStream inn=new DataInputStream(in);
				BufferedReader input= new BufferedReader(new InputStreamReader(inn));
				input.readLine();
				String s="";
				while((s=input.readLine())!=null){
					String[] t=s.split("	");
					if(!is_olfactory(t[0])){
						counts.add(Integer.parseInt(t[2]));
					}
				}
				input.close();
				
				
			}
			catch(Exception e){
				System.out.println(e);
			}
			ArrayList<int[]> histogram=histogram(counts);
			Collections.sort(histogram,comp);
			double[][] solution=minimize_parallel(50, histogram,  24, 100);//(150, histogram,  24, 100)
			
			
//			System.out.println(entity_name[k]);
//			for (int i=0;i<solution.length;i++){
//				System.out.println(ou(solution[i])+", "+(i+1));
//			}
//			for (int i=0;i<histogram.size();i++){
//				System.out.println(histogram.get(i)[0]+"	"+histogram.get(i)[1]+"	"+model1_density(histogram.get(i)[0],solution[0][0],solution[0][1])+"	"+model2_density(histogram.get(i)[0],solution[1][0],solution[1][1])+"	"+model3_density(histogram.get(i)[0],solution[2][0],solution[2][1],solution[2][2],solution[2][3])+"	"+model4_density(histogram.get(i)[0],solution[3][0],solution[3][1],solution[3][2],solution[3][3])+"	"+model5_density(histogram.get(i)[0],solution[4][0],solution[4][1],solution[4][2],solution[4][3],solution[4][4])+"	"+model6_density(histogram.get(i)[0],solution[5][0],solution[5][1],solution[5][2],solution[5][3],solution[5][4]));
//			}
			
			try{
				FileWriter out=new FileWriter(args[0]+"CBASE/Parameters_Summary/Parameters"+entity_name[k]+".txt");
				BufferedWriter output=new BufferedWriter(out);
				for (int i=0;i<solution.length;i++){
					output.write(ou(solution[i])+", "+(i+1));
					output.newLine();
				}
				output.close();
			}
			catch(Exception e){
				System.out.println(e);
			}
		}
	
	}
	
	public static int[] index_header(String[] header, String[] ideal_header){
		int[] indices=new int[ideal_header.length];
		for (int i=0;i<ideal_header.length;i++){
			int index=-1;
			for (int j=0;j<header.length;j++){
				if(header[j].equals(ideal_header[i])){
					index=j;
					break;
				}
			}
			indices[i]=index;
		}
		return indices;
	}
	
	public static boolean contains (String s, ArrayList<String> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).equals(s)){
				return true;
			}
		}
		return false;
	}
	
	
	public static double[] improve_short(int mod_C,  ArrayList<int[]> histogram, double[] start, double[][] bounds, int iterations){
		double[] a=minimize_neg_ln_L_short(discard_last(start),mod_C,bounds,histogram,iterations);
		double value=model_short(histogram,discard_last(a),mod_C);
		for (int i=0;i<100;i++){
			a=minimize_neg_ln_L_short(discard_last(a),mod_C,bounds,histogram,iterations);
			double value_new=model_short(histogram,discard_last(a),mod_C);
			if(value==value_new){
				break;
			}
			value=value_new;
			//System.out.println(ou(a));
		}
		return a;
	}
	
	public static double[] improve_long(int mod_C,  ArrayList<int[]> histogram, double[] start, double[][] bounds, int iterations){
		double[] a=minimize_neg_ln_L_long(discard_last(start),mod_C,bounds,histogram,iterations);
		double value=model_long(histogram,discard_last(a),mod_C);
		for (int i=0;i<100;i++){
			a=minimize_neg_ln_L_long(discard_last(a),mod_C,bounds,histogram,iterations);
			double value_new=model_long(histogram,discard_last(a),mod_C);
			if(value==value_new){
				break;
			}
			value=value_new;
			//System.out.println(ou(a));
		}
		return a;
	}
	
	
	public static double[] discard_last(double[] a){
		double[] b=new double[a.length-1];
		for (int i=0;i<a.length-1;i++){
			b[i]=a[i];
		}
		return b;
	}
	
	public static double [] minimize_neg_ln_L_short(double[] start, int mod, double[][] bounds, ArrayList<int[]> histo ,int iterations){//
		//bounds are currently ignored
		//double[][][] standard_bounds={{{0,100},{0,100}},{{0,100},{0,100}},{{0,100},{0,100},{0,100},{0,1}},{{0,100},{0,100},{0,100},{0,1}},{{0,100},{0,100},{0,100},{0,100},{0,1}},{{0,100},{0,100},{0,100},{0,100},{0,1}}};
		double[][][] standard_bounds={{{0.01,100},{0.01,100}},{{0.01,100},{0.01,100}},{{0.01,100},{0.01,100},{0.01,100},{0.01,1}},{{0.01,100},{0.01,100},{0.01,100},{0.01,1}},{{0.01,100},{0.01,100},{0.01,100},{0.01,100},{0.01,1}},{{0.01,100},{0.01,100},{0.01,100},{0.01,100},{0.01,1}}};
		
		
		//double[] x=start;
		//double[] g=gradient(x);
		ArrayList<double[]> g=new ArrayList<double[]>();//stores the last m+1 entries
		ArrayList<double[]> x=new ArrayList<double[]>();//stores the last m+1 entries
		
		x.add(start);
		g.add(gradient_short(x.get(x.size()-1),mod, histo));
//		System.out.println((0)+"	"+ou(x.get(x.size()-1))+"	"+model(histo,x.get(x.size()-1),mod));
		
		double value_best_solution=model_short(histo,start,mod);
		double[] best_solution=clone(start);
		
		outer:
		for (int iter=0;iter<iterations;iter++){
			double[] z=new double[0];
			//double scaling=1;
			if(x.size()>1){
				double[] q=clone(g.get(g.size()-1));
				double [] alpha=new double[x.size()-1];
				for (int i=1;i<x.size();i++){
					double[] y=diff(g.get(x.size()-i),g.get(x.size()-i-1));
					double[] s=diff(x.get(x.size()-i),x.get(x.size()-i-1));
					double rho=1/product_tv_v(y,s);
					alpha[x.size()-i-1]=rho*product_tv_v(s,q);
					q=diff(q,product(alpha[i-1],y));
				}
				
				
				double[] y=diff(g.get(g.size()-1),g.get(g.size()-2));
				double[] s=diff(x.get(x.size()-1),x.get(x.size()-2));
				
				//double scaling=product_tv_v(s,y)/sq_norm(y);
				double[][] h_k_0=product(1.0/sq_norm(y),product_v_tv(y,s));
				z=product(h_k_0,q);
				for (int i=1;i<x.size();i++){
					y=diff(g.get(i),g.get(i-1));
					s=diff(x.get(i),x.get(i-1));
					double rho=1/product_tv_v(y,s);
					double beta=rho*product_tv_v(y,z);
					z=add(z,product(alpha[i-1]-beta,s));
				}
				//System.out.println("z: "+ou(z));
				if(is_nan(z)){
					break outer;
				}
				
				z=product(1/Math.sqrt(sq_norm(z)),z);
				
				double min=Double.MAX_VALUE;
				double k_min=model_short(histo,x.get(x.size()-1),mod);
				for (double k=0.01;k<=1;k+=0.01){
					double m=model_short(histo,diff(x.get(x.size()-1),product(k,z)),mod);
					if(m<min&&in_bound(diff(x.get(x.size()-1),product(k,z)),bounds)){
						min=m;
						k_min=k;
					}
				}
				if(min!=Double.MAX_VALUE){
					z=product(k_min,z);
				}
				else{
					min=Double.MAX_VALUE;
					k_min=model_short(histo,x.get(x.size()-1),mod);
					for (double k=0.0001;k<=0.01;k+=0.0001){
						double m=model_short(histo,diff(x.get(x.size()-1),product(k,z)),mod);
						if(m<min&&in_bound(diff(x.get(x.size()-1),product(k,z)),bounds)){
							min=m;
							k_min=k;
						}
					}
					if(min!=Double.MAX_VALUE){
						z=product(k_min,z);
					}
					else{
						break outer;
					}
					
					
				}
			}
			else{
				if(mod==1||mod==2){
					z=new double[]{0.01,0.01};
				}
				else if(mod==3||mod==4){
					z= new double[]{0.01,0.01,0.01,0.01};
				}
				else if(mod==5||mod==6){
					z= new double[]{0.01,0.01,0.01,0.01,0.01};
				}
			}
			
		
			x.add(diff(x.get(x.size()-1),z));//add scaling factor
			g.add(gradient_short(x.get(x.size()-1), mod, histo));
//			System.out.println((iter+1)+"	"+ou(x.get(x.size()-1))+"	"+model(histo,x.get(x.size()-1),mod));
			
			double val=model_short(histo,x.get(x.size()-1),mod);
			if(val<value_best_solution&&val>0&&in_standard_bounds(x.get(x.size()-1),standard_bounds[mod-1])){
				value_best_solution=val;
				best_solution=clone(x.get(x.size()-1));
			}
			
		}
		
		
		
		return concat(best_solution,value_best_solution);//concat(x.get(x.size()-1),model(histo,x.get(x.size()-1),mod));
	}
	
	public static double [] minimize_neg_ln_L_long(double[] start, int mod, double[][] bounds, ArrayList<int[]> histo ,int iterations){//
		//bounds are currently ignored
		//double[][][] standard_bounds={{{0,100},{0,100}},{{0,100},{0,100}},{{0,100},{0,100},{0,100},{0,1}},{{0,100},{0,100},{0,100},{0,1}},{{0,100},{0,100},{0,100},{0,100},{0,1}},{{0,100},{0,100},{0,100},{0,100},{0,1}}};
		double[][][] standard_bounds={{{0.01,100},{0.01,100}},{{0.01,100},{0.01,100}},{{0.01,100},{0.01,100},{0.01,100},{0.01,1}},{{0.01,100},{0.01,100},{0.01,100},{0.01,1}},{{0.01,100},{0.01,100},{0.01,100},{0.01,100},{0.01,1}},{{0.01,100},{0.01,100},{0.01,100},{0.01,100},{0.01,1}}};
		
		
		//double[] x=start;
		//double[] g=gradient(x);
		ArrayList<double[]> g=new ArrayList<double[]>();//stores the last m+1 entries
		ArrayList<double[]> x=new ArrayList<double[]>();//stores the last m+1 entries
		
		x.add(start);
		g.add(gradient_long(x.get(x.size()-1),mod, histo));
//		System.out.println((0)+"	"+ou(x.get(x.size()-1))+"	"+model(histo,x.get(x.size()-1),mod));
		
		double value_best_solution=model_long(histo,start,mod);
		double[] best_solution=clone(start);
		
		outer:
		for (int iter=0;iter<iterations;iter++){
			double[] z=new double[0];
			//double scaling=1;
			if(x.size()>1){
				double[] q=clone(g.get(g.size()-1));
				double [] alpha=new double[x.size()-1];
				for (int i=1;i<x.size();i++){
					double[] y=diff(g.get(x.size()-i),g.get(x.size()-i-1));
					double[] s=diff(x.get(x.size()-i),x.get(x.size()-i-1));
					double rho=1/product_tv_v(y,s);
					alpha[x.size()-i-1]=rho*product_tv_v(s,q);
					q=diff(q,product(alpha[i-1],y));
				}
				
				
				double[] y=diff(g.get(g.size()-1),g.get(g.size()-2));
				double[] s=diff(x.get(x.size()-1),x.get(x.size()-2));
				
				//double scaling=product_tv_v(s,y)/sq_norm(y);
				double[][] h_k_0=product(1.0/sq_norm(y),product_v_tv(y,s));
				z=product(h_k_0,q);
				for (int i=1;i<x.size();i++){
					y=diff(g.get(i),g.get(i-1));
					s=diff(x.get(i),x.get(i-1));
					double rho=1/product_tv_v(y,s);
					double beta=rho*product_tv_v(y,z);
					z=add(z,product(alpha[i-1]-beta,s));
				}
				//System.out.println("z: "+ou(z));
				if(is_nan(z)){
					break outer;
				}
				
				z=product(1/Math.sqrt(sq_norm(z)),z);
				
				double min=Double.MAX_VALUE;
				double k_min=model_long(histo,x.get(x.size()-1),mod);
				for (double k=0.01;k<=1;k+=0.01){
					double m=model_long(histo,diff(x.get(x.size()-1),product(k,z)),mod);
					if(m<min&&in_bound(diff(x.get(x.size()-1),product(k,z)),bounds)){
						min=m;
						k_min=k;
					}
				}
				if(min!=Double.MAX_VALUE){
					z=product(k_min,z);
				}
				else{
					min=Double.MAX_VALUE;
					k_min=model_long(histo,x.get(x.size()-1),mod);
					for (double k=0.0001;k<=0.01;k+=0.0001){
						double m=model_long(histo,diff(x.get(x.size()-1),product(k,z)),mod);
						if(m<min&&in_bound(diff(x.get(x.size()-1),product(k,z)),bounds)){
							min=m;
							k_min=k;
						}
					}
					if(min!=Double.MAX_VALUE){
						z=product(k_min,z);
					}
					else{
						break outer;
					}
					
					
				}
			}
			else{
				if(mod==1||mod==2){
					z=new double[]{0.01,0.01};
				}
				else if(mod==3||mod==4){
					z= new double[]{0.01,0.01,0.01,0.01};
				}
				else if(mod==5||mod==6){
					z= new double[]{0.01,0.01,0.01,0.01,0.01};
				}
			}
			
		
			x.add(diff(x.get(x.size()-1),z));//add scaling factor
			g.add(gradient_long(x.get(x.size()-1), mod, histo));
//			System.out.println((iter+1)+"	"+ou(x.get(x.size()-1))+"	"+model(histo,x.get(x.size()-1),mod));
			
			double val=model_long(histo,x.get(x.size()-1),mod);
			if(val<value_best_solution&&val>0&&in_standard_bounds(x.get(x.size()-1),standard_bounds[mod-1])){
				value_best_solution=val;
				best_solution=clone(x.get(x.size()-1));
			}
			
		}
		
		
		
		return concat(best_solution,value_best_solution);//concat(x.get(x.size()-1),model(histo,x.get(x.size()-1),mod));
	}
	
	
	/*
	public static double [] minimize_neg_ln_L(double[] start, int mod, double[][] bounds, ArrayList<int[]> histo ,int iterations, double step){//
		//bounds are currently ignored
		
		//double[] x=start;
		//double[] g=gradient(x);
		ArrayList<double[]> g=new ArrayList<double[]>();//stores the last m+1 entries
		ArrayList<double[]> x=new ArrayList<double[]>();//stores the last m+1 entries
		
		x.add(start);
		g.add(gradient(x.get(x.size()-1),mod, histo));
//		System.out.println((0)+"	"+ou(x.get(x.size()-1))+"	"+model(histo,x.get(x.size()-1),mod));
		
		double value_best_solution=model(histo,start,mod);
		double[] best_solution=clone(start);
		
		outer:
		for (int iter=0;iter<iterations;iter++){
			double[] z=new double[0];
			//double scaling=1;
			if(x.size()>1){
				double[] q=clone(g.get(g.size()-1));
				double [] alpha=new double[x.size()-1];
				for (int i=1;i<x.size();i++){
					double[] y=diff(g.get(x.size()-i),g.get(x.size()-i-1));
					double[] s=diff(x.get(x.size()-i),x.get(x.size()-i-1));
					double rho=1/product_tv_v(y,s);
					alpha[x.size()-i-1]=rho*product_tv_v(s,q);
					q=diff(q,product(alpha[i-1],y));
				}
				
				
				double[] y=diff(g.get(g.size()-1),g.get(g.size()-2));
				double[] s=diff(x.get(x.size()-1),x.get(x.size()-2));
				
				//double scaling=product_tv_v(s,y)/sq_norm(y);
				double[][] h_k_0=product(1.0/sq_norm(y),product_v_tv(y,s));
				z=product(h_k_0,q);
				for (int i=1;i<x.size();i++){
					y=diff(g.get(i),g.get(i-1));
					s=diff(x.get(i),x.get(i-1));
					double rho=1/product_tv_v(y,s);
					double beta=rho*product_tv_v(y,z);
					z=add(z,product(alpha[i-1]-beta,s));
				}
				//System.out.println("z: "+ou(z));
				if(is_nan(z)){
					break outer;
				}
				
				z=product(1/Math.sqrt(sq_norm(z)),z);
				
				double min=Double.MAX_VALUE;
				double k_min=model(histo,x.get(x.size()-1),mod);
				for (double k=0.01*step;k<=1*step;k+=0.01*step){
					double m=model(histo,diff(x.get(x.size()-1),product(k,z)),mod);
					if(m<min&&in_bound(diff(x.get(x.size()-1),product(k,z)),bounds)){
						min=m;
						k_min=k;
					}
				}
				if(min!=Double.MAX_VALUE){
					z=product(k_min,z);
				}
				else{
					min=Double.MAX_VALUE;
					k_min=model(histo,x.get(x.size()-1),mod);
					for (double k=0.0001*step;k<=0.01*step;k+=0.0001*step){
						double m=model(histo,diff(x.get(x.size()-1),product(k,z)),mod);
						if(m<min&&in_bound(diff(x.get(x.size()-1),product(k,z)),bounds)){
							min=m;
							k_min=k;
						}
					}
					if(min!=Double.MAX_VALUE){
						z=product(k_min,z);
					}
					else{
						break outer;
					}
					
					
				}
			}
			else{
				if(mod==1||mod==2){
					z=new double[]{0.01,0.01};
				}
				else if(mod==3||mod==4){
					z= new double[]{0.01,0.01,0.01,0.01};
				}
				else if(mod==5||mod==6){
					z= new double[]{0.01,0.01,0.01,0.01,0.01};
				}
			}
			
		
			x.add(diff(x.get(x.size()-1),z));//add scaling factor
			g.add(gradient(x.get(x.size()-1), mod, histo));
//			System.out.println((iter+1)+"	"+ou(x.get(x.size()-1))+"	"+model(histo,x.get(x.size()-1),mod));
			
			double val=model(histo,x.get(x.size()-1),mod);
			if(val<value_best_solution){
				value_best_solution=val;
				best_solution=clone(x.get(x.size()-1));
			}
			
		}
		
		
		
		return concat(best_solution,value_best_solution);//concat(x.get(x.size()-1),model(histo,x.get(x.size()-1),mod));
	}
	*/
	
	public static boolean is_nan(double[] a){
		for (int i=0;i<a.length;i++){
			if(Double.isNaN(a[i])){
				return true;
			}
		}
		return false;
	}
	
	
	public static boolean in_standard_bounds(double [] x, double[][] bounds){
		for (int i=0;i<bounds.length;i++){
			if(x[i]<=bounds[i][0]||bounds[i][1]<x[i]){
				return false;
			}
		}
		return true;
	}
	
	
	
	public static boolean in_bound(double [] x, double[][] bounds){
		for (int i=0;i<x.length;i++){
			if(x[i]<bounds[i][0]||bounds[i][1]<x[i]){
				return false;
			}
		}
		return true;
	}
	
	public static String ou(double[] a){
		String s=""+a[0];
		for (int i=1;i<a.length;i++){
			s=s+", "+a[i];
		}
		return s;
	}
	
	public static double[] concat(double [] x, double y){
		double [] xx=new double[x.length+1];
		for (int i=0;i<x.length;i++){
			xx[i]=x[i];
		}
		xx[xx.length-1]=y;
		return xx;
	}
	
	
	public static double[] gradient_short(double[] v, int mode, ArrayList<int[]> histo){
		if(mode==1){
			return model1_short_diff(histo,v[0],v[1]);
		}
		else if(mode==2){
			return model2_short_diff(histo,v[0],v[1]);
		} 
		else if(mode==3){
			return model3_short_diff(histo,v[0],v[1],v[2],v[3]);
		} 
		else if(mode==4){
			return model4_short_diff(histo,v[0],v[1],v[2],v[3]);
		} 
		else if(mode==5){
			return model5_short_diff(histo,v[0],v[1],v[2],v[3],v[4]);
		}
		else if(mode==6){
			return model6_short_diff(histo,v[0],v[1],v[2],v[3],v[4]);
		} 
		else{
			return new double[0];
		}
	}
	
	
	public static double[] gradient_long(double[] v, int mode, ArrayList<int[]> histo){
		double delta=0.0001;
		if(mode==1){
			return gradient1(v,histo,delta);
			
		}
		else if(mode==2){
			return gradient2(v,histo,delta);
		} 
		else if(mode==3){
			return gradient3(v,histo,delta);
		} 
		else if(mode==4){
			return gradient4(v,histo,0.0000001);
		} 
		else if(mode==5){
			return gradient5(v,histo,0.0000001);
		}
		else if(mode==6){
			return gradient6(v,histo,delta);
		} 
		else{
			return new double[0];
		}
	}
	
	
	public static double[] gradient1(double[] v, ArrayList<int[]> histo, double delta){
		return new double[] {(model1(histo,v[0]+delta,v[1])-model1(histo,v[0]-delta,v[1]))/(2*delta),(model1(histo,v[0],v[1]+delta)-model1(histo,v[0],v[1]-delta))/(2*delta)};
	}
	
	public static double[] gradient2(double[] v, ArrayList<int[]> histo, double delta){
		return new double[] {(model2(histo,v[0]+delta,v[1])-model2(histo,v[0]-delta,v[1]))/(2*delta),(model2(histo,v[0],v[1]+delta)-model2(histo,v[0],v[1]-delta))/(2*delta)};
	}
	
	public static double[] gradient3(double[] v, ArrayList<int[]> histo, double delta){
		return new double[] {(model3(histo,v[0]+delta,v[1],v[2],v[3])-model3(histo,v[0]-delta,v[1],v[2],v[3]))/(2*delta),(model3(histo,v[0],v[1]+delta,v[2],v[3])-model3(histo,v[0],v[1]-delta,v[2],v[3]))/(2*delta),(model3(histo,v[0],v[1],v[2]+delta,v[3])-model3(histo,v[0],v[1],v[2]-delta,v[3]))/(2*delta),(model3(histo,v[0],v[1],v[2],v[3]+delta)-model3(histo,v[0],v[1],v[2],v[3]-delta))/(2*delta)};
	}
	
	public static double[] gradient4(double[] v, ArrayList<int[]> histo, double delta){
		return new double[] {(model4(histo,v[0]+delta,v[1],v[2],v[3])-model4(histo,v[0]-delta,v[1],v[2],v[3]))/(2*delta),(model4(histo,v[0],v[1]+delta,v[2],v[3])-model4(histo,v[0],v[1]-delta,v[2],v[3]))/(2*delta),(model4(histo,v[0],v[1],v[2]+delta,v[3])-model4(histo,v[0],v[1],v[2]-delta,v[3]))/(2*delta),(model4(histo,v[0],v[1],v[2],v[3]+delta)-model4(histo,v[0],v[1],v[2],v[3]-delta))/(2*delta)};
	}
	
	public static double[] gradient5(double[] v, ArrayList<int[]> histo, double delta){
		return new double[] {(model5(histo,v[0]+delta,v[1],v[2],v[3],v[4])-model5(histo,v[0]-delta,v[1],v[2],v[3],v[4]))/(2*delta),(model5(histo,v[0],v[1]+delta,v[2],v[3],v[4])-model5(histo,v[0],v[1]-delta,v[2],v[3],v[4]))/(2*delta),(model5(histo,v[0],v[1],v[2]+delta,v[3],v[4])-model5(histo,v[0],v[1],v[2]-delta,v[3],v[4]))/(2*delta),(model5(histo,v[0],v[1],v[2],v[3]+delta,v[4])-model5(histo,v[0],v[1],v[2],v[3]-delta,v[4]))/(2*delta),(model5(histo,v[0],v[1],v[2],v[3],v[4]+delta)-model5(histo,v[0],v[1],v[2],v[3],v[4]-delta))/(2*delta)};
	}
	
	public static double[] gradient6(double[] v, ArrayList<int[]> histo, double delta){
		return new double[] {(model6(histo,v[0]+delta,v[1],v[2],v[3],v[4])-model6(histo,v[0]-delta,v[1],v[2],v[3],v[4]))/(2*delta),(model6(histo,v[0],v[1]+delta,v[2],v[3],v[4])-model6(histo,v[0],v[1]-delta,v[2],v[3],v[4]))/(2*delta),(model6(histo,v[0],v[1],v[2]+delta,v[3],v[4])-model6(histo,v[0],v[1],v[2]-delta,v[3],v[4]))/(2*delta),(model6(histo,v[0],v[1],v[2],v[3]+delta,v[4])-model6(histo,v[0],v[1],v[2],v[3]-delta,v[4]))/(2*delta),(model6(histo,v[0],v[1],v[2],v[3],v[4]+delta)-model6(histo,v[0],v[1],v[2],v[3],v[4]-delta))/(2*delta)};
	}
	
	
	public static double sq_norm(double[] v){
		double sum=0;
		for (int i=0;i<v.length;i++){
			sum+=v[i]*v[i];
		}
		return sum;
	}
	
	public static double[] clone(double[] a){
		double[] b=new double[a.length];
		for (int i=0;i<a.length;i++){
			b[i]=a[i];
		}
		return b;
	}
	
	public static double[] diff(double[] v, double[] w){
		double[] diff=new double[v.length];
		for (int i=0;i<v.length;i++){
			diff[i]=v[i]-w[i];
		}
		return diff;
	}
	
	public static double[] add(double[] v, double[] w){
		double[] add=new double[v.length];
		for (int i=0;i<v.length;i++){
			add[i]=v[i]+w[i];
		}
		return add;
	}
	
	public static double[] product(double[][] m, double[] v){
		double[] w=new double[m.length];
		for (int i=0;i<m.length;i++){
			for (int j=0;j<m[i].length;j++){
				w[i]+=m[i][j]*v[j];
			}
		}
		return w;
	}
	
	public static double[][] product_v_tv(double[] v, double[] w){
		double [][] m=new double[v.length][w.length];
		for (int i=0;i<v.length;i++){
			for (int j=0;j<w.length;j++){
				m[i][j]=v[i]*w[j];
			}
		}
		return m;
	}
	
	public static double product_tv_v(double[] v, double[] w){
		double m=0;
		for (int i=0;i<v.length;i++){
			for (int j=0;j<w.length;j++){
				m+=v[i]*w[j];
			}
		}
		return m;
	}

	public static double[][] product(double a, double[][] m){
		double[][] n=new double[m.length][m[0].length];
		for (int i=0;i<m.length;i++){
			for (int j=0;j<m[i].length;j++){
				n[i][j]=a*m[i][j];
			}
		}
		return n;
	}
	
	public static double[] product(double a, double[] v){
		double[] w=new double[v.length];
		for (int i=0;i<v.length;i++){
			w[i]=a*v[i];
		}
		return w;
	}
	
	
	public static ArrayList<int[]> histogram(ArrayList<Integer> counts){
		ArrayList<int[]> histogram=new ArrayList<int[]>();
		for(int i=0;i<counts.size();i++){
			int ii=index(counts.get(i).intValue(),histogram);
			if(ii==-1){
				histogram.add(new int[]{counts.get(i),1});
			}
			else{
				histogram.get(ii)[1]++;
			}
		}
		return histogram;
	}
	
	public static int index(int a, ArrayList<int[]> b){
		for (int i=0;i<b.size();i++){
			if(b.get(i)[0]==a){
				return i;
			}
		}
		return -1;
	}
	
	public static boolean is_olfactory(String s){
		
		return s.length()>=3&&s.substring(0,2).equals("OR")&&isInteger(s.substring(2,3));
	}
	
	public static boolean isInteger(String s){
		try{
			int a=Integer.parseInt(s);
			return true;
		}
		catch(Exception e){
			return false;
		}
	}
	
	
	//public static double[] minimize_neg_ln_L(double[] start, int mod_C, double[][] bounds){
	//	return new double[0];
	//}
	
	public static double[][] minimize_parallel(int rep_no, ArrayList<int[]> histo, int no_cpu, int iterations){
		ArrayList<double[]>[] solutions=new ArrayList[6];
		for (int mod=1;mod<=6;mod++){
			//if(mod<6){
			//	continue;
			//}
			solutions[mod-1]=minimize_parallel(mod,rep_no,histo,no_cpu,iterations);
		}
		/*
		for (int j=0;j<solutions.length;j++){
			if(j<5){
				continue;
			}
			for (int k=0;k<solutions[j].size();k++){
				for (int l=0;l<solutions[j].get(k).length;l++){
					System.out.print("	"+solutions[j].get(k)[l]);
				}
				System.out.println();
			}
		}
		System.out.println(solutions[5].size());
		System.exit(0);*/
		
		double[][][] standard_bounds={{{0.01,100},{0.01,100}},{{0.01,100},{0.01,100}},{{0.01,100},{0.01,100},{0.01,100},{0.01,1}},{{0.01,100},{0.01,100},{0.01,100},{0.01,1}},{{0.01,100},{0.01,100},{0.01,100},{0.01,100},{0.01,1}},{{0.01,100},{0.01,100},{0.01,100},{0.01,100},{0.01,1}}};
		
		/*
		SubthreadImprove[][] subthreads=new SubthreadImprove[solutions.length][];
		for (int i=0;i<solutions.length;i++){
			subthreads[i]=new SubthreadImprove[solutions[i].size()];
			for (int j=0;j<solutions[i].size();j++){
				subthreads[i][j]=new SubthreadImprove(i+1,  histo, solutions[i].get(j), standard_bounds[i],  iterations);
			}
		}
		execute(subthreads,no_cpu);
		*/
		double[][][] subthreads=new double[solutions.length][][];
		for (int i=0;i<solutions.length;i++){
			subthreads[i]=new double[solutions[i].size()][];
			for (int j=0;j<solutions[i].size();j++){
				subthreads[i][j]=run_improve(i+1,  histo, solutions[i].get(j), standard_bounds[i],  iterations);
			}
		}
		
		
		double[][] solutions2=new double[6][];
		for (int i=0;i<subthreads.length;i++){
			double min=Double.MAX_VALUE;
			int j_min=-1;
			for (int j=0;j<subthreads[i].length;j++){
				if(subthreads[i][j][subthreads[i][j].length-1]<min){
					min=subthreads[i][j][subthreads[i][j].length-1];
					j_min=j;
				}
			}
			if(j_min!=-1){
				solutions2[i]=subthreads[i][j_min];	
			}
			else{
				solutions2[i]=new double[0];
			}
			
		}
		/*
		for (int i=0;i<subthreads.length;i++){
			double min=Double.MAX_VALUE;
			int j_min=-1;
			for (int j=0;j<subthreads[i].length;j++){
				if(subthreads[i][j].solution[subthreads[i][j].solution.length-1]<min){
					min=subthreads[i][j].solution[subthreads[i][j].solution.length-1];
					j_min=j;
				}
			}
			if(j_min!=-1){
				solutions2[i]=subthreads[i][j_min].solution;	
			}
			else{
				solutions2[i]=new double[0];
			}
			
		}*/
		return solutions2;
	}
	
	
	public static ArrayList<double[]> minimize_parallel(int mod_C, int rep_no, ArrayList<int[]> histo, int no_cpu, int iterations){
		//double[][][] standard_bounds={{{0,100},{0,100}},{{0,100},{0,100}},{{0,100},{0,100},{0,100},{0,1}},{{0,100},{0,100},{0,100},{0,1}},{{0,100},{0,100},{0,100},{0,100},{0,1}},{{0,100},{0,100},{0,100},{0,100},{0,1}}};
		double[][][] standard_bounds={{{0.01,100},{0.01,100}},{{0.01,100},{0.01,100}},{{0.01,100},{0.01,100},{0.01,100},{0.01,1}},{{0.01,100},{0.01,100},{0.01,100},{0.01,1}},{{0.01,100},{0.01,100},{0.01,100},{0.01,100},{0.01,1}},{{0.01,100},{0.01,100},{0.01,100},{0.01,100},{0.01,1}}};
		
		
		
		ArrayList<double[]> solutions=new ArrayList<double[]>();
		if(mod_C==1){
			//double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] low_b=new double[]{Math.pow(10,-2)*uniform(1,3),Math.pow(10,-2)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2)};
			
			//cur_min_res=new double[]{0,0,Double.MAX_VALUE};
			
			
			double[][] threads=new double[rep_no][];
			for (int i=0;i<rep_no;i++){
				threads[i] =  run(new double[]{uniform(0.02,10), uniform(0.02,10)}, 1, new double[][]{{low_b[0],up_b[0]}, {low_b[1],up_b[1]}}, histo,iterations);	
			}
			for (int i=0;i<threads.length;i++){
				if(in_standard_bounds(threads[i],standard_bounds[0])){
					//System.out.println((i+1)+"	"+ou(p_res));
					if (threads[i][2]>0 && threads[i][2]<Double.MAX_VALUE){
						solutions.add(new double[]{threads[i][0],threads[i][1],threads[i][2]});	
					}
				}
			}
			
			
		}
		else if(mod_C==2){
			//double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] low_b=new double[]{Math.pow(10,-2)*uniform(1,3),Math.pow(10,-2)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2)};
			//cur_min_res=new double[]{0,0,Double.MAX_VALUE};
			
			double[][] threads=new double[rep_no][];
			for (int i=0;i<rep_no;i++){
				threads[i]=run(new double[]{uniform(0.02,10), uniform(0.02,10)},2, new double[][]{{low_b[0],up_b[0]}, {low_b[1],up_b[1]}}, histo,iterations);
			}
			for (int i=0;i<threads.length;i++){
				if(in_standard_bounds(threads[i],standard_bounds[1])){
					//System.out.println((i+1)+"	"+ou(p_res));
					if (threads[i][2]>0 && threads[i][2]<Double.MAX_VALUE){
						solutions.add(new double[]{threads[i][0],threads[i][1],threads[i][2]});	
					}
				}
			}
			
		}
		else if(mod_C==3){
			//double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] low_b=new double[]{Math.pow(10,-2)*uniform(1,3),Math.pow(10,-2)*uniform(1,3),Math.pow(10,-2)*uniform(1,3),Math.pow(10,-2)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2)};
			//cur_min_res=new double[]{0,0,0,0,Double.MAX_VALUE};
			
			double[][] threads=new double[rep_no][];
			for (int i=0;i<rep_no;i++){
				//threads[i]=new Subthread(new double[]{uniform(0.02,10), uniform(0.02,10), uniform(0.02,10), uniform(2*Math.pow(10, -5),0.95)},  3, new double[][]{{low_b[0],up_b[0]},{low_b[1],up_b[1]},{low_b[2],up_b[2]},{low_b[3],0.9999}},histo,iterations);
				threads[i]=run(new double[]{uniform(0.02,10), uniform(0.02,10), uniform(0.02,10), uniform(2*Math.pow(10, -2),0.95)},  3, new double[][]{{low_b[0],up_b[0]},{low_b[1],up_b[1]},{low_b[2],up_b[2]},{low_b[3],0.9999}},histo,iterations);
				
			}
			for (int i=0;i<threads.length;i++){
				if(in_standard_bounds(threads[i],standard_bounds[2])){
					if (threads[i][4]>0 && threads[i][4]<Double.MAX_VALUE){
						solutions.add(new double[]{threads[i][0],threads[i][1],threads[i][2],threads[i][3],threads[i][4]});
					}	
				}	
			}
			
		}
		else if(mod_C==4){
			//double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] low_b=new double[]{Math.pow(10,-2)*uniform(1,3),Math.pow(10,-2)*uniform(1,3),Math.pow(10,-2)*uniform(1,3),Math.pow(10,-2)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2)};
			//cur_min_res=new double[]{0,0,0,0,Double.MAX_VALUE};
			
			double[][] threads=new double[rep_no][];
			for (int i=0;i<rep_no;i++){
				//threads[i]=new Subthread(new double[]{uniform(0.02,10), uniform(0.02,10), uniform(0.02,10), uniform(2*Math.pow(10,-5),0.95)}, 4, new double[][]{{low_b[0],up_b[0]},{low_b[1],up_b[1]},{low_b[2],up_b[2]},{low_b[3],0.9999}},histo,iterations);
				threads[i]=run(new double[]{uniform(0.02,10), uniform(0.02,10), uniform(0.02,10), uniform(2*Math.pow(10,-2),0.95)}, 4, new double[][]{{low_b[0],up_b[0]},{low_b[1],up_b[1]},{low_b[2],up_b[2]},{low_b[3],0.9999}},histo,iterations);
				
			}
			for (int i=0;i<threads.length;i++){
				if(in_standard_bounds(threads[i],standard_bounds[3])){
					if (threads[i][4]>0 && threads[i][4]<Double.MAX_VALUE){
						solutions.add(new double[]{threads[i][0],threads[i][1],threads[i][2],threads[i][3],threads[i][4]});
					}	
				}
			}
			
		}
		else if(mod_C==5){
			//double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] low_b=new double[]{Math.pow(10,-2)*uniform(1,3),Math.pow(10,-2)*uniform(1,3),Math.pow(10,-2)*uniform(1,3),Math.pow(10,-2)*uniform(1,3),Math.pow(10,-2)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2)};
			//cur_min_res=new double[]{0,0,0,0,0,Double.MAX_VALUE};
			
			double[][] threads=new double[rep_no][];
			for (int i=0;i<rep_no;i++){
				//threads[i]=new Subthread(new double[]{uniform(0.02,10), uniform(0.02,5), uniform(0.02,10), uniform(0.02,10.),uniform(2*Math.pow(10, -5),0.95)},  5, new double[][]{{low_b[0],up_b[0]},{low_b[1],up_b[1]},{low_b[2],up_b[2]},{low_b[3],up_b[3]},{low_b[4],0.9999}},histo,iterations);
				threads[i]=run(new double[]{uniform(0.02,10), uniform(0.02,5), uniform(0.02,10), uniform(0.02,10.),uniform(2*Math.pow(10, -2),0.95)},  5, new double[][]{{low_b[0],up_b[0]},{low_b[1],up_b[1]},{low_b[2],up_b[2]},{low_b[3],up_b[3]},{low_b[4],0.9999}},histo,iterations);
				
			}
			for (int i=0;i<threads.length;i++){
				if(in_standard_bounds(threads[i],standard_bounds[4])){
					if (threads[i][5]>0 && threads[i][5]<Double.MAX_VALUE){
						solutions.add(new double[]{threads[i][0],threads[i][1],threads[i][2],threads[i][3],threads[i][4],threads[i][5]});
					}
				}
			}
			
		}	
		else if(mod_C==6){
			//double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] low_b=new double[]{Math.pow(10,-2)*uniform(1,3),Math.pow(10,-2)*uniform(1,3),Math.pow(10,-2)*uniform(1,3),Math.pow(10,-2)*uniform(1,3),Math.pow(10,-2)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2)};
			//cur_min_res=new double[]{0,0,0,0,0,Double.MAX_VALUE};
			
			double[][] threads=new double[rep_no][];
			for (int i=0;i<rep_no;i++){
				//threads[i]=new Subthread(new double[]{uniform(0.02,10), uniform(0.02,5), uniform(0.02,10), uniform(0.02,10.),uniform(2*Math.pow(10, -5),0.95)},  6, new double[][]{{low_b[0],up_b[0]},{low_b[1],up_b[1]},{low_b[2],up_b[2]},{low_b[3],up_b[3]},{low_b[4],0.9999}},histo,iterations);
				threads[i]=run(new double[]{uniform(0.02,10), uniform(0.02,5), uniform(0.02,10), uniform(0.02,10.),uniform(2*Math.pow(10, -2),0.95)},  6, new double[][]{{low_b[0],up_b[0]},{low_b[1],up_b[1]},{low_b[2],up_b[2]},{low_b[3],up_b[3]},{low_b[4],0.9999}},histo,iterations);
				
			}
			for (int i=0;i<threads.length;i++){
				//System.out.println(threads[i][0]+"	"+threads[i][1]+"	"+threads[i][2]+"	"+threads[i][3]+"	"+threads[i][4]+"	"+threads[i][5]);
				//System.out.println(model_long(histo,threads[i],mod_C)+"	"+model_short(histo,threads[i],mod_C));
				if(in_standard_bounds(threads[i],standard_bounds[5])){
					if (threads[i][5]>0 && threads[i][5]<Double.MAX_VALUE){
						solutions.add(new double[]{threads[i][0],threads[i][1],threads[i][2],threads[i][3],threads[i][4],threads[i][5]});
					}
				}
			}
		}
		Comparator<double[]> comp=(double[] a, double[] b)->{
			return new Double(a[a.length-1]).compareTo(b[b.length-1]);
		};
		Collections.sort(solutions,comp);
		
		ArrayList<double[]> solutions_final=new ArrayList<double[]>();
		for (int i=0;i<Math.min(5, solutions.size());i++){
			solutions_final.add(solutions.get(i));
		}	
		return solutions_final;
	}
	
	/*
	public static void execute(Subthread[] threads, int no_cpu){

		int no_undone=threads.length;
		int no_running=0;
		while(no_undone>0){
			no_undone=0;
			no_running=0;
			for (int i=0;i<threads.length;i++){
				if(!threads[i].done){
					no_undone++;
				}
				if(threads[i].running){
					no_running++;
				}
			}
			if(no_undone==0){
				break;
			}
			if(no_running<no_cpu){
				for (int i=0;i<threads.length;i++){
					if(!threads[i].done&&!threads[i].running){
						threads[i].start();
						no_running++;
					}
					if(no_running>=no_cpu){
						break;
					}
				}
			}
			try{
				Thread.sleep(500);
			}
			catch(Exception e){
				System.out.println(e);
			}
		}
	}*/
	
	/*
	public static void execute(SubthreadImprove[][] threads, int no_cpu){

		int no_undone=threads.length;
		int no_running=0;
		while(no_undone>0){
			no_undone=0;
			no_running=0;
			for (int i=0;i<threads.length;i++){
				for (int j=0;j<threads[i].length;j++){
					if(!threads[i][j].done){
						no_undone++;
					}
					if(threads[i][j].running){
						no_running++;
					}
				}
				
			}
			if(no_undone==0){
				break;
			}
			if(no_running<no_cpu){
				for (int i=0;i<threads.length;i++){
					for (int j=0;j<threads[i].length;j++){
						if(!threads[i][j].done&&!threads[i][j].running){
							threads[i][j].start();
							no_running++;
						}
						if(no_running>=no_cpu){
							break;
						}
					}
				}
			}
			try{
				Thread.sleep(500);
			}
			catch(Exception e){
				System.out.println(e);
			}
		}
	}*/
	/*
	private static class SubthreadImprove extends Thread{
		volatile double[] solution=new double[0];
		volatile boolean done=false;
		volatile boolean running=false;
		int mod_C=-1;
		ArrayList<int[]> histogram=new ArrayList<int[]>();
		double[] start=new double[0];
		double[][] bounds=new double[0][0];
		int iterations=0;
		
		public SubthreadImprove(int mod_C,  ArrayList<int[]> histogram, double[] start, double[][] bounds, int iterations){
			this.mod_C=mod_C;
			this.histogram=histogram;
			this.start=start;
			this.bounds=bounds;
			this.iterations=iterations;
		}
		
		public void run(){
			running=true;
			done=false;
			solution=improve_short(mod_C, histogram, start, bounds, iterations);//TODO: short vs long
			running=false;
			done=true;
		}
	}*/
	
	public static double[]  run_improve(int mod_C, ArrayList<int[]> histogram,double[]  start,double[][] bounds, int iterations){
		return improve_short(mod_C, histogram, start, bounds, iterations);//TODO: short vs long
	}
	
	/*
	private static class Subthread extends Thread{
		volatile double[] solution=new double[0];
		volatile boolean done=false;
		volatile boolean running=false;
		 start=new double[0];
		 mod=-1;
		bounds=new double[0][0];
		 histo=new ArrayList<int[]>();
		 iterations=0;
		public Subthread(double[] start, int mod, double[][] bounds, ArrayList<int[]> histo, int iterations){
			this.start=start;
			this.mod=mod;
			this.bounds=bounds;
			this.histo=histo;
			this.iterations=iterations;
		}
		public void run(){
			running=true;
			done=false;
			solution=minimize_neg_ln_L_short( start,  mod, bounds,  histo, iterations);//TODO: short vs long
			running=false;
			done=true;
		}
	}*/
	
	public static double[] run(double[] start, int mod,double[][]  bounds, ArrayList<int[]> histo, int iterations){
		return minimize_neg_ln_L_short(start,  mod, bounds,  histo, iterations);//TODO: short vs long
	}
	
	
	/*
	public static double[] minimizeee(int mod_C, int rep_no, ArrayList<int[]> histo, int iterations){
		double[] cur_min_res=new double[0];
		if(mod_C==1){
			double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2)};
			double[][] standard_bounds={{0,100},{0,100}};
			cur_min_res=new double[]{0,0,Double.MAX_VALUE};
			for (int i=0;i<rep_no;i++){
				double[] p_res = minimize_neg_ln_L(new double[]{uniform(0.02,10), uniform(0.02,10)}, 1, new double[][]{{low_b[0],up_b[0]}, {low_b[1],up_b[1]}}, histo,iterations);
				if(in_standard_bounds(p_res,standard_bounds)){
					System.out.println((i+1)+"	"+ou(p_res));
					if (p_res[2]>0 && p_res[2]<cur_min_res[2]){
						cur_min_res =new double[]{p_res[0],p_res[1],p_res[2]};	
					}
				}		
			}
			
		}
		else if(mod_C==2){
			double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2)};
			cur_min_res=new double[]{0,0,Double.MAX_VALUE};
			double[][] standard_bounds={{0,100},{0,100}};
			for (int i=0;i<rep_no;i++){
				double[] p_res = minimize_neg_ln_L(new double[]{uniform(0.02,10), uniform(0.02,10)},2, new double[][]{{low_b[0],up_b[0]}, {low_b[1],up_b[1]}}, histo,iterations);
				if(in_standard_bounds(p_res,standard_bounds)){
					System.out.println((i+1)+"	"+ou(p_res));
					if (p_res[2]>0 && p_res[2]<cur_min_res[2]){
						cur_min_res =new double[]{p_res[0],p_res[1],p_res[2]};	
					}
				}
				
			}
			
		}
		else if(mod_C==3){
			double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2)};
			cur_min_res=new double[]{0,0,0,0,Double.MAX_VALUE};
			double[][] standard_bounds={{0,100},{0,100},{0,100},{0,1}};
			for (int i=0;i<rep_no;i++){
				double[] p_res = minimize_neg_ln_L(new double[]{uniform(0.02,10), uniform(0.02,10), uniform(0.02,10), uniform(2*Math.pow(10, -5),0.95)},  3, new double[][]{{low_b[0],up_b[0]},{low_b[1],up_b[1]},{low_b[2],up_b[2]},{low_b[3],0.9999}},histo,iterations);
				if(in_standard_bounds(p_res,standard_bounds)){
					System.out.println((i+1)+"	"+ou(p_res));
					if (p_res[4]>0 && p_res[4]<cur_min_res[4]){
						cur_min_res =new double[]{p_res[0],p_res[1],p_res[2],p_res[3],p_res[4]};
					}
				}
				
			}
		}
		else if(mod_C==4){
			double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2)};
			cur_min_res=new double[]{0,0,0,0,Double.MAX_VALUE};
			double[][] standard_bounds={{0,100},{0,100},{0,100},{0,1}};
			for (int i=0;i<rep_no;i++){
				double[] p_res = minimize_neg_ln_L(new double[]{uniform(0.02,10), uniform(0.02,10), uniform(0.02,10), uniform(2*Math.pow(10,-5),0.95)}, 4, new double[][]{{low_b[0],up_b[0]},{low_b[1],up_b[1]},{low_b[2],up_b[2]},{low_b[3],0.9999}},histo,iterations);
				if(in_standard_bounds(p_res,standard_bounds)){
					System.out.println((i+1)+"	"+ou(p_res));
					if (p_res[4]>0 && p_res[4]<cur_min_res[4]){
						cur_min_res =new double[]{p_res[0],p_res[1],p_res[2],p_res[3],p_res[4]};
					}
				}
				
			}
		}
		else if(mod_C==5){
			double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2)};
			cur_min_res=new double[]{0,0,0,0,0,Double.MAX_VALUE};
			double[][] standard_bounds={{0,100},{0,100},{0,100},{0,100},{0,1}};
			for (int i=0;i<rep_no;i++){
				double[] p_res = minimize_neg_ln_L(new double[]{uniform(0.02,10), uniform(0.02,5), uniform(0.02,10), uniform(0.02,10.),uniform(2*Math.pow(10, -5),0.95)},  5, new double[][]{{low_b[0],up_b[0]},{low_b[1],up_b[1]},{low_b[2],up_b[2]},{low_b[3],up_b[3]},{low_b[4],0.9999}},histo,iterations);
				if(in_standard_bounds(p_res,standard_bounds)){
					System.out.println((i+1)+"	"+ou(p_res));
					if (p_res[5]>0 && p_res[5]<cur_min_res[5]){
						cur_min_res =new double[]{p_res[0],p_res[1],p_res[2],p_res[3],p_res[4],p_res[5]};
					}
				}
				
			}
		}	
		else if(mod_C==6){
			double[] low_b=new double[]{Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3),Math.pow(10,-5)*uniform(1,3)};
			double[] up_b=new double[]{50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2),50*uniform(1,2)};
			cur_min_res=new double[]{0,0,0,0,0,Double.MAX_VALUE};
			double[][] standard_bounds={{0,100},{0,100},{0,100},{0,100},{0,1}};
			for (int i=0;i<rep_no;i++){
				double[] p_res = minimize_neg_ln_L(new double[]{uniform(0.02,10), uniform(0.02,5), uniform(0.02,10), uniform(0.02,10.),uniform(2*Math.pow(10, -5),0.95)},  6, new double[][]{{low_b[0],up_b[0]},{low_b[1],up_b[1]},{low_b[2],up_b[2]},{low_b[3],up_b[3]},{low_b[4],0.9999}},histo,iterations);
				if(in_standard_bounds(p_res,standard_bounds)){
					System.out.println((i+1)+"	"+ou(p_res));
					if (p_res[5]>0 && p_res[5]<cur_min_res[5]){
						cur_min_res =new double[]{p_res[0],p_res[1],p_res[2],p_res[3],p_res[4],p_res[5]};
						
					}
				}
				
			}
		}
		return cur_min_res;
	}
	*/
	
	
	public static double model_long(ArrayList<int[]> histo, double[] param, int mod){
		if(mod==1){
			return model1(histo,param[0],param[1]);
		}
		else if(mod==2){
			return model2(histo,param[0],param[1]);
		}
		else if(mod==3){
			return model3(histo,param[0],param[1],param[2],param[3]);
		}
		else if(mod==4){
			return model4(histo,param[0],param[1],param[2],param[3]);
		}
		else if(mod==5){
			return model5(histo,param[0],param[1],param[2],param[3],param[4]);
		}
		else if(mod==6){
			return model6(histo,param[0],param[1],param[2],param[3],param[4]);
		}
		else{
			return 0;
		}
	}
	
	public static double model_short(ArrayList<int[]> histo, double[] param, int mod){
		if(mod==1){
			return model1_short(histo,param[0],param[1]);
		}
		else if(mod==2){
			return model2_short(histo,param[0],param[1]);
		}
		else if(mod==3){
			return model3_short(histo,param[0],param[1],param[2],param[3]);
		}
		else if(mod==4){
			return model4_short(histo,param[0],param[1],param[2],param[3]);
		}
		else if(mod==5){
			return model5_short(histo,param[0],param[1],param[2],param[3],param[4]);
		}
		else if(mod==6){
			return model6_short(histo,param[0],param[1],param[2],param[3],param[4]);
		}
		else{
			return 0;
		}
	}
	
	
	public static double model1(ArrayList<int[]> h, double a, double b){
		double sum=0;
		for (int i=0;i<h.size();i++){
			double x= (h.get(i)[0]*Math.log(b) + (-h.get(i)[0]-a)*Math.log(1 + b) + Gamma.logGamma(h.get(i)[0] + a) -  Gamma.logGamma(h.get(i)[0]+1) -  Gamma.logGamma(a));
			if(Double.isFinite(x)){
				sum+=x*h.get(i)[1];
			}
		}
		return -sum;
	}
	
	public static double model1_short(ArrayList<int[]> h, double a, double b){
		double sum=0;
		for (int i=0;i<h.size();i++){
			double x= (h.get(i)[0]*log(b) + (-h.get(i)[0]-a)*log(1 + b) + logGamma(h.get(i)[0] + a) - logGamma(h.get(i)[0]+1) -  logGamma(a));
			if(Double.isFinite(x)){
				sum+=x*h.get(i)[1];
			}
			
		}
		return -sum;
	}
	
	
	public static double[] model1_short_diff(ArrayList<int[]> h, double a, double b){
		double[] sum=new double[2];
		for (int i=0;i<h.size();i++){
			double x_1= ( -log(1 + b) + delta_logGamma(h.get(i)[0] + a)  -  delta_logGamma(a));
			double x_2= (h.get(i)[0]*delta_log(b) + (-h.get(i)[0]-a)*delta_log(1 + b) );
			
			if(Double.isFinite(x_1)&&Double.isFinite(x_2)){
				sum[0]-=x_1*h.get(i)[1];
				sum[1]-=x_2*h.get(i)[1];
			}
		}

		return sum;
	}
	
	
	
	
	public static double log(double x){
		
		if(Double.isNaN(x)||x<=1.1||(int)(x*1000)>=log_precompute.length){
			return Math.log(x);
		}
		else {
			return log_precompute[(int)(x*1000)];
		}
		
	}
	public static double logGamma(double x){
		if(Double.isNaN(x)||x<=1||(int)(x*1000)>=gamma_precompute.length){
			return Gamma.logGamma(x);
		}
		else{
			return gamma_precompute[(int)(x*1000)];
		}
	}
	public static double delta_logGamma(double x){
		if(Double.isNaN(x)||x<=1||(int)(x*1000)>=gamma_precompute.length){
			return (Gamma.logGamma(x+0.0001)-Gamma.logGamma(x-0.0001))/(2*0.0001);
		}
		else{
			return gamma_delta_precompute[(int)(x*1000)];
		}
		
	}
	public static double delta_log(double x){
		return 1.0/x;
	}
	public static double exp(double x){
		
		if(!Double.isNaN(x)&&0<=(int)((x+1000)*1000)&&(int)((x+1000)*1000)<exp_precompute.length){
			return exp_precompute[(int)((x+1000)*1000)];
		}
		else{
			return Math.exp(x);
		}
		//return Math.exp(x);
	}
	public static double besselk_diff1_short(double n, double x){
		/*
		if(Double.isNaN(x)||Double.isNaN(n)){
			return Double.NaN;
		}
		else if(0<=(int)((n+20)*100)&&(int)((n+20)*100)<kv_precompute_diff1[0].length&&0<=(int)(x*100)&&(int)(x*100)<kv_precompute_diff1.length){
			//System.out.println(kv_precompute_diff1[(int)(x*100)][(int)((n+20)*100)]+"	"+((Bessel.k(x, n+0.0001, false)-Bessel.k(x, n-0.0001, false))/(2*0.0001)));
			return kv_precompute_diff1[(int)(x*100)][(int)((n+20)*100)];
		}
		else{
			return (Bessel.k(x, n+0.0001, false)-Bessel.k(x, n-0.0001, false))/(2*0.0001);
		}*/
		/*
		int x1=(int)(x*100);
		int x2=(int)(x*100)+1;
		int n1=(int)((n+20)*100);
		int n2=(int)((n+20)*100)+1;
		
		if(Double.isNaN(x)||Double.isNaN(n)){
			return Double.NaN;
		}
		else if(0<=n1&&n2<kv_precompute_diff1[0].length&&0<=x1&&x2<kv_precompute_diff1.length){
			double d_x=(x-(double)(x1)/100.0)/((double)(x2)/100.0-(double)(x1)/100.0);
			double d_n=(n-((double)(n1)/100.0-20))/((double)(n2)/100.0-(double)(n1)/100.0);
			double a=Math.exp(Math.log(kv_precompute_diff1[x1][n1])*(1-d_x)*(1-d_n)+Math.log(kv_precompute_diff1[x1][n2])*(1-d_x)*(d_n)+Math.log(kv_precompute_diff1[x2][n1])*(d_x)*(1-d_n)+Math.log(kv_precompute_diff1[x2][n2])*d_x*d_n);
			//System.out.println(a+"	"+(Bessel.k(x, n+0.0001, false)-Bessel.k(x, n-0.0001, false))/(2*0.0001));
			if(Double.isNaN(a)||a==0){
				return (Bessel.k(x, n+0.0001, false)-Bessel.k(x, n-0.0001, false))/(2*0.0001);
				
			}
			else{
				return a;
			}
			
			//Bessel.k(x, n, false);
		}
		else{
			return (Bessel.k(x, n+0.0001, false)-Bessel.k(x, n-0.0001, false))/(2*0.0001);
		}*/		
		return (Bessel.k(x, n+0.0001, false)-Bessel.k(x, n-0.0001, false))/(2*0.0001);
		
	}
	public static double besselk_diff2_short(double n, double x){
		/*
		if(Double.isNaN(x)||Double.isNaN(n)){
			return Double.NaN;
		}
		else if(0<=(int)((n+20)*100)&&(int)((n+20)*100)<kv_precompute_diff2[0].length&&0<=(int)(x*100)&&(int)(x*100)<kv_precompute_diff2.length){
			return kv_precompute_diff2[(int)(x*100)][(int)((n+20)*100)];
		}
		else{
			return (Bessel.k(x+0.0001, n, false)-Bessel.k(x-0.0001, n, false))/(2*0.0001);
		}*/
		
		/*
		int x1=(int)(x*100);
		int x2=(int)(x*100)+1;
		int n1=(int)((n+20)*100);
		int n2=(int)((n+20)*100)+1;
		
		if(Double.isNaN(x)||Double.isNaN(n)){
			return Double.NaN;
		}
		else if(0<=n1&&n2<kv_precompute_diff2[0].length&&0<=x1&&x2<kv_precompute_diff2.length){
			double d_x=(x-(double)(x1)/100.0)/((double)(x2)/100.0-(double)(x1)/100.0);
			double d_n=(n-((double)(n1)/100.0-20))/((double)(n2)/100.0-(double)(n1)/100.0);
			double a=-Math.exp(Math.log(-kv_precompute_diff2[x1][n1])*(1-d_x)*(1-d_n)+Math.log(-kv_precompute_diff2[x1][n2])*(1-d_x)*(d_n)+Math.log(-kv_precompute_diff2[x2][n1])*(d_x)*(1-d_n)+Math.log(-kv_precompute_diff2[x2][n2])*d_x*d_n);
			
			if(Double.isNaN(a)||a==0){
				return (Bessel.k(x+0.0001, n, false)-Bessel.k(x-0.0001, n, false))/(2*0.0001);
				
			}
			else{
				return a;
			}
			//System.out.println(a+"	"+(Bessel.k(x+0.0001, n, false)-Bessel.k(x-0.0001, n, false))/(2*0.0001));
			
			//System.out.println(a+"	"+(Bessel.k(x, n+0.0001, false)-Bessel.k(x, n-0.0001, false))/(2*0.0001));
			//Bessel.k(x, n, false);
		}
		else{
			return (Bessel.k(x+0.0001, n, false)-Bessel.k(x-0.0001, n, false))/(2*0.0001);
		}*/
		return (Bessel.k(x+0.0001, n, false)-Bessel.k(x-0.0001, n, false))/(2*0.0001);
		
	}
	
	
	public static double besselk_short(double n, double x){
		/*
		if(Double.isNaN(x)||Double.isNaN(n)){
			return Double.NaN;
		}
		else if(0<=(int)((n+20)*100)&&(int)((n+20)*100)<kv_precompute[0].length&&0<=(int)(x*100)&&(int)(x*100)<kv_precompute.length){
			//System.out.println(n+"	"+x+"	"+kv_precompute[(int)(x*100)][(int)((n+20)*100)]+"	"+Bessel.k(x, n, false));
			return kv_precompute[(int)(x*100)][(int)((n+20)*100)];
		}
		else{
			return Bessel.k(x, n, false);
		}*/
		
		/*
		int x1=(int)(x*100);
		int x2=(int)(x*100)+1;
		int n1=(int)((n+20)*100);
		int n2=(int)((n+20)*100)+1;
		
		if(Double.isNaN(x)||Double.isNaN(n)){
			return Double.NaN;
		}
		else if(0<=n1&&n2<kv_precompute[0].length&&0<=x1&&x2<kv_precompute.length){
			double d_x=(x-(double)(x1)/100.0)/((double)(x2)/100.0-(double)(x1)/100.0);
			double d_n=(n-((double)(n1)/100.0-20))/((double)(n2)/100.0-(double)(n1)/100.0);
			double a=Math.exp(Math.log(kv_precompute[x1][n1])*(1-d_x)*(1-d_n)+Math.log(kv_precompute[x1][n2])*(1-d_x)*(d_n)+Math.log(kv_precompute[x2][n1])*(d_x)*(1-d_n)+Math.log(kv_precompute[x2][n2])*d_x*d_n);
			//System.out.println(a+"	"+Bessel.k(x, n, false));
			return a;//Bessel.k(x, n, false);
		}
		else{
			return Bessel.k(x, n, false);
		}*/
		return Bessel.k(x, n, false);
	}
	
	
	public static double besselk_long(double n, double x){
		return Bessel.k(x, n, false);
	}
	public static double delta_logGamma_long(double x){
		return (Gamma.logGamma(x+0.0001)-Gamma.logGamma(x-0.0001))/(2*0.0001);
	}
	
	
	
	
	public static double model2(ArrayList<int[]> h, double a, double b){
		double sum=0;
		for (int i=0;i<h.size();i++){
			double proxy=besselk_long(-h.get(i)[0] + a, 2*Math.sqrt(b));
			//System.out.println(h.get(i)[0]+"	"+proxy);
			if(proxy>0&&Double.isFinite(proxy)){
				double x= Math.log(2) + ((h.get(i)[0] + a)/2.0)*Math.log(b) + Math.log(proxy) - Gamma.logGamma(h.get(i)[0]+1) - Gamma.logGamma(a);
				if(Double.isFinite(x)){
					sum+=x*h.get(i)[1];
				}
				
			}
		}
		return -sum;		
	}
	
	public static double model2_short(ArrayList<int[]> h, double a, double b){
		double sum=0;
		for (int i=0;i<h.size();i++){
			double proxy=besselk_short(-h.get(i)[0] + a, 2*Math.sqrt(b));
			if(proxy>0&&Double.isFinite(proxy)){
				double x= log(2) + ((h.get(i)[0] + a)/2.0)*log(b) + log(proxy) - logGamma(h.get(i)[0]+1) - logGamma(a);
				if(Double.isFinite(x)){
					sum+=x*h.get(i)[1];
				}
				
			}
		}
		return -sum;			
	}
	
	public static double[] model2_short_diff(ArrayList<int[]> h, double a, double b){
		
		double[] sum=new double[2];
		for (int i=0;i<h.size();i++){
			double proxy=besselk_short(-h.get(i)[0] + a, 2*Math.sqrt(b));
			double proxy_diff_1=besselk_diff1_short(-h.get(i)[0] + a, 2*Math.sqrt(b));
			double proxy_diff_2=besselk_diff2_short(-h.get(i)[0] + a, 2*Math.sqrt(b));
			
			if(proxy>0&&Double.isFinite(proxy)){
				double x_1=  + 0.5*log(b) + delta_log(proxy)*proxy_diff_1  - delta_logGamma(a);
				double x_2= ((h.get(i)[0] + a)/2.0)*delta_log(b) + delta_log(proxy)*proxy_diff_2/Math.sqrt(b);
				
				if(Double.isFinite(x_1)&&Double.isFinite(x_2)){
					sum[0]-=x_1*h.get(i)[1];
					sum[1]-=x_2*h.get(i)[1];
				}
			}
		}
		
		return sum;			
		
	}
	
	
	
	public static double model3_short(ArrayList<int[]> h, double a, double b, double t, double w){
		double sum=0;
		for (int i=0;i<h.size();i++){
			double x=  log(exp( log(w * t) + (-1-h.get(i)[0])*log(1 + t) ) + exp( log(1-w) + h.get(i)[0]*log(b) + (-h.get(i)[0]-a)*log(1 + b) + logGamma(h.get(i)[0] + a) - logGamma(h.get(i)[0]+1) -logGamma(a) )) ;
			if(Double.isFinite(x)){
				sum+=x*h.get(i)[1];
			}
		}
		return -sum;
	}
	
	public static double model3(ArrayList<int[]> h, double a, double b, double t, double w){
		double sum=0;
		for (int i=0;i<h.size();i++){
			double x= Math.log( Math.exp( Math.log(w * t) + (-1-h.get(i)[0])*Math.log(1 + t) ) + Math.exp( Math.log(1-w) + h.get(i)[0]*Math.log(b) + (-h.get(i)[0]-a)*Math.log(1 + b) + Gamma.logGamma(h.get(i)[0] + a) - Gamma.logGamma(h.get(i)[0]+1) - Gamma.logGamma(a) ) );
			
			//System.out.println(h.get(i)[0]+"	"+x);
			if(Double.isFinite(x)){
				sum+=x*h.get(i)[1];
			}
		}
		return -sum;
	}
	
	public static double[] model3_short_diff(ArrayList<int[]> h, double a, double b, double t, double w){
		
		double[] sum=new double[4];
		
		for (int i=0;i<h.size();i++){
			double x1= delta_log( exp( log(w * t) + (-1-h.get(i)[0])*log(1 + t) ) + exp( log(1-w) + h.get(i)[0]*log(b) + (-h.get(i)[0]-a)*log(1 + b) + logGamma(h.get(i)[0] + a) - logGamma(h.get(i)[0]+1) -logGamma(a) ) );
			double x2=  exp( log(1-w) + h.get(i)[0]*log(b) + (-h.get(i)[0]-a)*log(1 + b) + logGamma(h.get(i)[0] + a) - logGamma(h.get(i)[0]+1) -logGamma(a) ) ;
			
			double x2_3=  exp( log(w * t) + (-1-h.get(i)[0])*log(1 + t) );
			double x2_4= exp(  (-1-h.get(i)[0])*log(1 + t) )*t  -exp(  h.get(i)[0]*log(b) + (-h.get(i)[0]-a)*log(1 + b) + logGamma(h.get(i)[0] + a) - logGamma(h.get(i)[0]+1) -logGamma(a) );
			
			double x3_1=  -log(1 + b) + delta_logGamma(h.get(i)[0] + a) -delta_logGamma(a);
			double x3_2= h.get(i)[0]*delta_log(b) + (-h.get(i)[0]-a)*delta_log(1 + b) ;
			double x3_3= delta_log(t) + (-1-h.get(i)[0])*delta_log(1 + t);
			
			
			double x_1=x1*x2*x3_1;
			double x_2=x1*x2*x3_2;
			double x_3=x1*x2_3*x3_3;
			double x_4=x1*x2_4;
			
			if(Double.isFinite(x_1)&&Double.isFinite(x_2)&&Double.isFinite(x_3)&&Double.isFinite(x_4)){
				sum[0]-=x_1*h.get(i)[1];
				sum[1]-=x_2*h.get(i)[1];
				sum[2]-=x_3*h.get(i)[1];
				sum[3]-=x_4*h.get(i)[1];
			}
			
		}

		return sum;
		
	}
	

	
	public static double model4_short(ArrayList<int[]> h, double a, double b, double t, double w){
		double sum=0;
		for (int i=0;i<h.size();i++){
				double x=  log(exp( log(w * t) + (-1 - h.get(i)[0])*log(1 + t) ) + exp( log(1-w) + log(2) + ((h.get(i)[0] + a)/2.0)*log(b) + log(besselk_short(-h.get(i)[0] + a, 2*Math.sqrt(b))) - logGamma(h.get(i)[0] + 1) - logGamma(a) ) );	
				if(Double.isFinite(x)){
					sum+=x*h.get(i)[1];
				}
			
			
		}
		return -sum;
	}

	public static double model4(ArrayList<int[]> h, double a, double b, double t, double w){
		double sum=0;
		for (int i=0;i<h.size();i++){
			double x= Math.log( Math.exp( Math.log(w * t) + (-1 - h.get(i)[0])*Math.log(1 + t) ) + Math.exp( Math.log(1-w) + Math.log(2) + ((h.get(i)[0] + a)/2.0)*Math.log(b) + Math.log(besselk_long(-h.get(i)[0] + a, 2*Math.sqrt(b))) - Gamma.logGamma(h.get(i)[0] + 1) - Gamma.logGamma(a) ) );	
			if(Double.isFinite(x)){
				sum+=x*h.get(i)[1];
			}
		}
		return -sum;
	}
	
	
	
	public static double model5(ArrayList<int[]> h, double a, double b, double g, double d, double w){
		double sum=0;
		for (int i=0;i<h.size();i++){
			double x= Math.log( Math.exp( Math.log(w) + h.get(i)[0]*Math.log(b) + (-h.get(i)[0]-a)*Math.log(1 + b) + Gamma.logGamma(h.get(i)[0] + a) - Gamma.logGamma(h.get(i)[0]+1) - Gamma.logGamma(a) ) + Math.exp( Math.log(1-w) + h.get(i)[0]*Math.log(d) + (-h.get(i)[0]-g)*Math.log(1 + d) + Gamma.logGamma(h.get(i)[0] + g) - Gamma.logGamma(h.get(i)[0]+1) - Gamma.logGamma(g) ) );
			if(Double.isFinite(x)){
				sum+=x*h.get(i)[1];
			}
		}
		return -sum;
	}
	
	public static double model5_short(ArrayList<int[]> h, double a, double b, double g, double d, double w){
		double sum=0;
		for (int i=0;i<h.size();i++){
			double x=  log(exp( log(w) + h.get(i)[0]*log(b) + (-h.get(i)[0]-a)*log(1 + b) + logGamma(h.get(i)[0] + a) - logGamma(h.get(i)[0]+1) - logGamma(a) ) + exp( log(1-w) + h.get(i)[0]*log(d) + (-h.get(i)[0]-g)*log(1 + d) + logGamma(h.get(i)[0] + g) - logGamma(h.get(i)[0]+1) - logGamma(g) )) ;
			if(Double.isFinite(x)){
				sum+=x*h.get(i)[1];
			}
		}
		return -sum;
	}
	
	public static double[] model4_short_diff(ArrayList<int[]> h, double a, double b, double t, double w){
		
		double[] sum=new double[4];
		for (int i=0;i<h.size();i++){
			double x1= delta_log( exp( log(w * t) + (-1 - h.get(i)[0])*log(1 + t) ) + exp( log(1-w) + log(2) + ((h.get(i)[0] + a)/2.0)*log(b) + log(besselk_short(-h.get(i)[0] + a, 2*Math.sqrt(b))) - logGamma(h.get(i)[0] + 1) - logGamma(a) ) );	
			double x2_3= exp( log(w * t) + (-1 - h.get(i)[0])*log(1 + t) );	
			double x2_4= exp( (-1 - h.get(i)[0])*log(1 + t) )*t  + exp(  log(2) + ((h.get(i)[0] + a)/2.0)*log(b) + log(besselk_short(-h.get(i)[0] + a, 2*Math.sqrt(b))) - logGamma(h.get(i)[0] + 1) - logGamma(a) )*(-1) ;	
			double x2= exp( log(1-w) + log(2) + ((h.get(i)[0] + a)/2.0)*log(b) + log(besselk_short(-h.get(i)[0] + a, 2*Math.sqrt(b))) - logGamma(h.get(i)[0] + 1) - logGamma(a) ) ;
				
			double x3_1= 0.5*log(b) + delta_log(besselk_short(-h.get(i)[0] + a, 2*Math.sqrt(b)))*besselk_diff1_short(-h.get(i)[0] + a, 2*Math.sqrt(b)) - delta_logGamma(a)  ;	
			double x3_2= ((h.get(i)[0] + a)/2.0)*delta_log(b) + delta_log(besselk_short(-h.get(i)[0] + a, 2*Math.sqrt(b)))*besselk_diff2_short(-h.get(i)[0] + a, 2*Math.sqrt(b))/Math.sqrt(b)   ;	
			double x3_3= delta_log( t) + (-1 - h.get(i)[0])*delta_log(1 + t);	
			
			double x_1=x1*x2*x3_1;
			double x_2=x1*x2*x3_2;
			double x_3=x1*x2_3*x3_3;
			double x_4=x1*x2_4;
			
			
			
			if(Double.isFinite(x_1)&&Double.isFinite(x_2)&&Double.isFinite(x_3)&&Double.isFinite(x_4)){
				sum[0]-=x_1*h.get(i)[1];
				sum[1]-=x_2*h.get(i)[1];
				sum[2]-=x_3*h.get(i)[1];
				sum[3]-=x_4*h.get(i)[1];
			}
			
		}
		return sum;
	}
	
	public static double[] model5_short_diff(ArrayList<int[]> h, double a, double b, double g, double d, double w){
		double[] sum=new double[5];
		for (int i=0;i<h.size();i++){
			double x1  = delta_log( exp( log(w) + h.get(i)[0]*log(b) + (-h.get(i)[0]-a)*log(1 + b) + logGamma(h.get(i)[0] + a) - logGamma(h.get(i)[0]+1) - logGamma(a) ) + exp( log(1-w) + h.get(i)[0]*log(d) + (-h.get(i)[0]-g)*log(1 + d) + logGamma(h.get(i)[0] + g) - logGamma(h.get(i)[0]+1) - logGamma(g) ) );
			double x2_12= exp( log(w) + h.get(i)[0]*log(b) + (-h.get(i)[0]-a)*log(1 + b) + logGamma(h.get(i)[0] + a) - logGamma(h.get(i)[0]+1) - logGamma(a) );
			double x2_34= exp( log(1-w) + h.get(i)[0]*log(d) + (-h.get(i)[0]-g)*log(1 + d) + logGamma(h.get(i)[0] + g) - logGamma(h.get(i)[0]+1) - logGamma(g) ) ;
			double x2_5= exp(  h.get(i)[0]*log(b) + (-h.get(i)[0]-a)*log(1 + b) + logGamma(h.get(i)[0] + a) - logGamma(h.get(i)[0]+1) - logGamma(a) ) -exp(  h.get(i)[0]*log(d) + (-h.get(i)[0]-g)*log(1 + d) + logGamma(h.get(i)[0] + g) - logGamma(h.get(i)[0]+1) - logGamma(g) );
			
			
			double x3_1= -log(1 + b) + delta_logGamma(h.get(i)[0] + a)  - delta_logGamma(a) ;
			double x3_2= h.get(i)[0]*delta_log(b) + (-h.get(i)[0]-a)*delta_log(1 + b)  ;
			double x3_3= -log(1 + d) + delta_logGamma(h.get(i)[0] + g) - delta_logGamma(g)  ;
			double x3_4= h.get(i)[0]*delta_log(d) + (-h.get(i)[0]-g)*delta_log(1 + d) ;
			
			
			double x_1=x1*x2_12*x3_1;
			double x_2=x1*x2_12*x3_2;
			double x_3=x1*x2_34*x3_3;
			double x_4=x1*x2_34*x3_4;
			double x_5=x1*x2_5;
			
			if(Double.isFinite(x_1)&&Double.isFinite(x_2)&&Double.isFinite(x_3)&&Double.isFinite(x_4)&&Double.isFinite(x_5)){
				sum[0]-=x_1*h.get(i)[1];
				sum[1]-=x_2*h.get(i)[1];
				sum[2]-=x_3*h.get(i)[1];
				sum[3]-=x_4*h.get(i)[1];
				sum[4]-=x_5*h.get(i)[1];
			}
			
		}
		
		
		return sum;
		
		
	}
	
	
	public static double model6(ArrayList<int[]> h, double a, double b, double g, double d, double w){
		double sum=0;
		for (int i=0;i<h.size();i++){
			double x= Math.log( (w * Math.exp(h.get(i)[0]*Math.log(b) + (-h.get(i)[0]-a)*Math.log(1 + b) + Gamma.logGamma(h.get(i)[0] + a) - Gamma.logGamma(h.get(i)[0]+1) - Gamma.logGamma(a))) + ((1-w) * Math.exp( Math.log(2) + ((h.get(i)[0] + g)/2.0)*Math.log(d) + Math.log(besselk_long(-h.get(i)[0] + g, 2*Math.sqrt(d))) - Gamma.logGamma(h.get(i)[0] + 1) - Gamma.logGamma(g) ) ));		
			
			if(Double.isFinite(x)){
				sum+=x*h.get(i)[1];
			}
				
			
		}
		return -sum;
	}
	
	public static double model6_short(ArrayList<int[]> h, double a, double b, double g, double d, double w){
		double sum=0;
		for (int i=0;i<h.size();i++){
			//if (h.get(i)[0]>25){
				double x=  (w * exp(h.get(i)[0]*log(b) + (-h.get(i)[0]-a)*log(1 + b) + logGamma(h.get(i)[0] + a) - logGamma(h.get(i)[0]+1) - logGamma(a))) + ((1-w) * exp( log(2) + ((h.get(i)[0] + g)/2.0)*log(d) + log(besselk_short(-h.get(i)[0] + g, 2*Math.sqrt(d))) - logGamma(h.get(i)[0] + 1) - logGamma(g) ) );		
				if(Double.isFinite(x)){
					sum+=log(x)*h.get(i)[1];
				}
				
			
		}
		return -sum;
	}
	
	public static double[] model6_short_diff(ArrayList<int[]> h, double a, double b, double g, double d, double w){
		double[] sum=new double[5];
		
		for (int i=0;i<h.size();i++){
			double x1  = delta_log( (w * exp(h.get(i)[0]*log(b) + (-h.get(i)[0]-a)*log(1 + b) + logGamma(h.get(i)[0] + a) - logGamma(h.get(i)[0]+1) - logGamma(a))) + ((1-w) * exp( log(2) + ((h.get(i)[0] + g)/2.0)*log(d) + log(besselk_short(-h.get(i)[0] + g, 2*Math.sqrt(d))) - logGamma(h.get(i)[0] + 1) - logGamma(g) ) ));		
			
			double x2_12= w * exp(h.get(i)[0]*log(b) + (-h.get(i)[0]-a)*log(1 + b) + logGamma(h.get(i)[0] + a) - logGamma(h.get(i)[0]+1) - logGamma(a)) ;		
			double x2_34= (1-w) * exp( log(2) + ((h.get(i)[0] + g)/2.0)*log(d) + log(besselk_short(-h.get(i)[0] + g, 2*Math.sqrt(d))) - logGamma(h.get(i)[0] + 1) - logGamma(g)) ;		
			double x2_5= exp(h.get(i)[0]*log(b) + (-h.get(i)[0]-a)*log(1 + b) + logGamma(h.get(i)[0] + a) - logGamma(h.get(i)[0]+1) - logGamma(a)) -exp( log(2) + ((h.get(i)[0] + g)/2.0)*log(d) + log(besselk_short(-h.get(i)[0] + g, 2*Math.sqrt(d))) - logGamma(h.get(i)[0] + 1) - logGamma(g)) ;		
			
			double x3_1= -log(1 + b) + delta_logGamma(h.get(i)[0] + a) - delta_logGamma(a) ;		
			double x3_2= h.get(i)[0]*delta_log(b) + (-h.get(i)[0]-a)*delta_log(1 + b) ;		
			double x3_3=  0.5*log(d) + delta_log(besselk_short(-h.get(i)[0] + g, 2*Math.sqrt(d)))*besselk_diff1_short(-h.get(i)[0] + g, 2*Math.sqrt(d))  - delta_logGamma(g) ;		
			double x3_4=  ((h.get(i)[0] + g)/2.0)*delta_log(d) + delta_log(besselk_short(-h.get(i)[0] + g, 2*Math.sqrt(d)))*besselk_diff2_short(-h.get(i)[0] + g, 2*Math.sqrt(d))/Math.sqrt(d) ;		
			
			
			double x_1=x1*x2_12*x3_1;
			double x_2=x1*x2_12*x3_2;
			double x_3=x1*x2_34*x3_3;
			double x_4=x1*x2_34*x3_4;
			double x_5=x1*x2_5;
			
			
			if(Double.isFinite(x_1)&&Double.isFinite(x_2)&&Double.isFinite(x_3)&&Double.isFinite(x_4)&&Double.isFinite(x_5)){
				sum[0]-=x_1*h.get(i)[1];
				sum[1]-=x_2*h.get(i)[1];
				sum[2]-=x_3*h.get(i)[1];
				sum[3]-=x_4*h.get(i)[1];
				sum[4]-=x_5*h.get(i)[1];
				
				
			}


		}
		
		
		return sum;
		
	}
	
	
	public static double uniform (double a, double b){
		return a+Math.random()*(b-a);
	}
	
	public static double model1_density(int h, double a, double b){
		 return Math.exp(h*Math.log(b) + (-h-a)*Math.log(1 + b) + Gamma.logGamma(h + a) -  Gamma.logGamma(h+1) -  Gamma.logGamma(a));
	}
	
	public static double model2_density(int h, double a, double b){//TODO: bessel
		double proxy=besselk_long(-h + a, 2*Math.sqrt(b));
		return Math.exp(Math.log(2) + ((h + a)/2.0)*Math.log(b) + Math.log(proxy) - Gamma.logGamma(h+1) - Gamma.logGamma(a));
	}
	
	public static double model3_density(int h, double a, double b, double t, double w){
		return Math.exp(Math.log( Math.exp( Math.log(w * t) + (-1-h)*Math.log(1 + t) ) + Math.exp( Math.log(1-w) + h*Math.log(b) + (-h-a)*Math.log(1 + b) + Gamma.logGamma(h + a) - Gamma.logGamma(h+1) - Gamma.logGamma(a) ) ));
	}
	
	public static double model4_density(int h, double a, double b, double t, double w){
		return Math.exp(Math.log( Math.exp( Math.log(w * t) + (-1 - h)*Math.log(1 + t) ) + Math.exp( Math.log(1-w) + Math.log(2) + ((h + a)/2.0)*Math.log(b) + Math.log(besselk_long(-h + a, 2*Math.sqrt(b))) - Gamma.logGamma(h + 1) - Gamma.logGamma(a) ) ));	
	}
	
	public static double model5_density(int h, double a, double b, double g, double d, double w){
		return  Math.exp(Math.log( Math.exp( Math.log(w) + h*Math.log(b) + (-h-a)*Math.log(1 + b) + Gamma.logGamma(h + a) - Gamma.logGamma(h+1) - Gamma.logGamma(a) ) + Math.exp( Math.log(1-w) + h*Math.log(d) + (-h-g)*Math.log(1 + d) + Gamma.logGamma(h + g) - Gamma.logGamma(h+1) - Gamma.logGamma(g) ) ));
	}
	
	public static double model6_density(int h, double a, double b, double g, double d, double w){
		return Math.exp(Math.log( (w * Math.exp(h*Math.log(b) + (-h-a)*Math.log(1 + b) + Gamma.logGamma(h + a) - Gamma.logGamma(h+1) - Gamma.logGamma(a))) + ((1-w) * Math.exp( Math.log(2) + ((h + g)/2.0)*Math.log(d) + Math.log(besselk_long(-h + g, 2*Math.sqrt(d))) - Gamma.logGamma(h + 1) - Gamma.logGamma(g) ) )));		
	}
	
	/*
	public static double kv(double n, double x){
		return Bessel.k(x, n, false);
	}*/
}
