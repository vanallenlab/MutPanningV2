/************************************************************           
 * MutPanning 												*
 * 															*   
 * Author:		Felix Dietlein								*   
 *															*   
 * Copyright:	(C) 2019 									*   
 *															*   
 * License:		BSD-3-Clause open source license			*  
 *															*   
 * Summary: This script clusters the count vectors 			*
 *			separately for each cancer entity. 				*
 * 															*   
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
//import java.io.File;

import org.apache.commons.math3.special.Gamma;

public class ClusteringEntity {

	static String[] entities=new String[0];
	static String file_type="";
	static String file_affinity="";
	static String file_samples="";
	static String file_reference="";
	static String file_out="";
	
	//static double[][] distance_type=new double[0][0];
	//static double[][] distance_context=new double[0][0];
	//static ArrayList<String> names_type=new ArrayList<String>();
	//static ArrayList<int[]> type=new ArrayList<int[]>();
	//static ArrayList<int[][][]> affinity=new ArrayList<int[][][]>();
	//static ArrayList<ArrayList<Integer>> cluster_context=new ArrayList<ArrayList<Integer>>();
	//static ArrayList<ArrayList<Integer>> cluster_type=new ArrayList<ArrayList<Integer>>();
	
	//static double[] vector_context=new double[0];
	//static double[] vector_type=new double[0];
	static double[][][] ref=new double[20][2][4];
	static String[] index_header_samples={"ID","Sample","Cohort"};
	
	
	
	static double[] gamma=new double[50000+1];
	
	/*
	 * argument0: root file
	 * argument1: sample annotaiton file
	 * argument2: no. processors to distribute workload
	 */
	
	public static void main(String[] args){
		
		file_type=args[0]+"AffinityCounts/TypeCount.txt";//"C:\\Users\\Administrator/Dropbox/AffinityCounts/CosmicCount.txt";
		file_affinity=args[0]+"AffinityCounts/AffinityCount.txt";//"C:\\Users\\Administrator/Dropbox/AffinityCounts/AffinityCount.txt";
		file_samples=args[1];
		file_reference=args[2]+"FileReferenceCount.txt";
		file_out=args[0]+"Clustering/";
		
		if(!new File(file_out).exists()){
			new File(file_out).mkdir();
		}
		
		
		//Determine all entity names, as clustering is performed separately for each Step this is needed to coordinate the order
		try{
			FileInputStream in=new FileInputStream(file_samples);
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
			entities=new String[aa.size()];
			for (int i=0;i<aa.size();i++){
				entities[i]=aa.get(i);
			}
		}
		catch(Exception e){
			System.out.println(e);
		}
		
		
	//	static String[] entities={"AdenoidCystic","Bladder","Blood","Brain","Breast","Cervix","Cholangio","Colorectal","Endometrium","Gastroesophageal","HeadNeck","KidneyClear","KidneyNonClear","Liver","LungAD","LungSC","LungSCLC","Lymph","Ovarian","Pancreas","Pheochromocytoma","Pleura","Prostate","Sarcoma","Skin","TesticularGermCell","Thymus","Thyroid","UvealMelanoma"};

		
		
		
		//save the Gamma values as they are needed frequently later, to save some time
		for (int i=0;i<gamma.length;i++){
			gamma[i]=Gamma.logGamma((double)(i)/10.0);
		}
		
		
		try{
			
			//counts in the reference genome,needed as prior in the distance metrics
			String s="";
			FileInputStream in=new FileInputStream(file_reference);
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			input.readLine();
			for (int i=0;i<2;i++){
				for (int j=0;j<4;j++){
					String[] t=input.readLine().split("	");
					for (int k=0;k<20;k++){
						if(k<10){
							ref[k][i][j]=Integer.parseInt(t[k+2]);
						}
						else {
							ref[k][i][j]=Integer.parseInt(t[k+3]);
						}
					}
				}
			}
			input.close();
			
			//normalization step reference 
			for (int i=0;i<ref.length;i++){
				for (int j=0;j<ref[i].length;j++){
					double sum=0;
					for (int k=0;k<ref[i][j].length;k++){
						sum+=ref[i][j][k];
					}
					for (int k=0;k<ref[i][j].length;k++){
						ref[i][j][k]/=sum;
					}
				}
			}
			
			/*
			for (int i=0;i<ref.length;i++){
				for (int j=0;j<ref[i].length;j++){
					System.out.println(ref[i][j][0]+"	"+ref[i][j][1]+"	"+ref[i][j][2]+"	"+ref[i][j][3]);
				}
			}
			System.exit(0);
			*/
			
			//for every entity clustering is first performed separately. clusters are merged
			//between entities later. as running time is n*n this safes a lot of time  
			for (int aa=0;aa<entities.length;aa++){
				System.out.println(entities[aa]);
				//reset all the variables
				ArrayList<String> names_type=new ArrayList<String>();
				ArrayList<int[]> type=new ArrayList<int[]>();
				ArrayList<int[][][]> affinity=new ArrayList<int[][][]>();
				
				//vector_context=new double[0];
				//vector_type=new double[0];
				
				
				
				//extract all samples needed in this step
				ArrayList<String> selected=new ArrayList<String>();
				in=new FileInputStream(file_samples);
				inn=new DataInputStream(in);
				input= new BufferedReader(new InputStreamReader(inn));
				int[] index_header=index_header(input.readLine().split("	"),index_header_samples);
				while((s=input.readLine())!=null){
					if(s.split("	")[index_header[2]].equals(entities[aa])){
						selected.add(s.split("	")[index_header[1]]);
					}
				}
				input.close();
				
				in=new FileInputStream(file_type);
				inn=new DataInputStream(in);
				input= new BufferedReader(new InputStreamReader(inn));
				input.readLine();
				
				FileInputStream in2=new FileInputStream(file_affinity);
				DataInputStream inn2=new DataInputStream(in2);
				BufferedReader input2= new BufferedReader(new InputStreamReader(inn2));
				input2.readLine();
				
				//read count vectors of the particular entity which will be clustered
				while((s=input.readLine())!=null){
					String[] t=s.split("	");
					String[] t2=input2.readLine().split("	");
					if(!t[0].equals(t2[0])){
						int a=1/0;
					}
					if(contains(t[0],selected)){
						
						int[] c=new int[6];
						for (int i=0;i<c.length;i++){
							c[i]=Integer.parseInt(t[i+1]);
						}
						int[][][] c2=new int[20][6][4];
						int n=0;
						for (int i=0;i<c2.length;i++){
							for (int j=0;j<c2[i].length;j++){
								for (int k=0;k<c2[i][j].length;k++){
									c2[i][j][k]=Integer.parseInt(t2[n+1]);
									n++;
								}
							}
						}
						
						if(sum(c)>0){
							names_type.add(t[0]);
							type.add(c);
							affinity.add(c2);
						}
						
					}
				}
				input.close();
				input2.close();
				
				
				double[][] distance_type=initialize_type(type);//new double[0][0];
				double[][] distance_context=initialize_context(affinity);//new double[0][0];
				
				/*
				for (int i=0;i<distance_type.length;i++){
					for (int j=0;j<distance_type[i].length;j++){
						System.out.print("	"+distance_type[i][j]);
					}
					System.out.println();
				}*/
				//System.exit(0);
				
				/*
				for (int i=0;i<distance_context.length;i++){
					for (int j=0;j<distance_context[i].length;j++){
						System.out.print("	"+distance_context[i][j]);
					}
					System.out.println();
				}
				System.exit(0);
				*/
				
				//distance_type=new double[type.size()][type.size()];
				//distance_context=new double[affinity.size()][affinity.size()];
				
				//compute the initial nxn distance metrics. this is a very compute intensive step
				// so that this computation step is performed in parallel 
				
				
				//initialize();
			

				//transform static matrix into a flexible ArrayList, which allows more flexible
				// rearrangements and deletions throughout clustering
				ArrayList<ArrayList<Double>> dist_type=new ArrayList<ArrayList<Double>>();
				for (int i=0;i<distance_type.length;i++){
					dist_type.add(new ArrayList<Double>());
					for (int j=0;j<distance_type[i].length;j++){
						dist_type.get(i).add(distance_type[i][j]);
					}
				}

				
				ArrayList<ArrayList<Double>> dist_context=new ArrayList<ArrayList<Double>>();
				for (int i=0;i<distance_context.length;i++){
					dist_context.add(new ArrayList<Double>());
					for (int j=0;j<distance_context[i].length;j++){
						dist_context.get(i).add(distance_context[i][j]);
					}
				}

				
				ArrayList<ArrayList<Integer>> cluster_type=new ArrayList<ArrayList<Integer>>();
				ArrayList<int[]> cluster_type_sum=new ArrayList<int[]>();
				for (int i=0;i<type.size();i++){
					ArrayList<Integer> c=new ArrayList<Integer>();
					c.add(i);
					cluster_type.add(c);
					cluster_type_sum.add(clone(type.get(i)));
				}
				
				
				//a very straighforward hierarchical clustering. always find the step with the minimal distance
				//and merge the clusters and update the distance metrics in each step
				//notice that the max is set to 0 (i.e. the ratio =1). hence, this procedure terminates when
				//the ratio falls below 1.
				//note that the update is performed on 30CPUs in parallel to speed up the process
				while(cluster_type.size()>0){
					int i_max=-1;
					int j_max=-1;
					double max=0;
					for (int i=0;i<cluster_type.size();i++){
						for (int j=i+1;j<cluster_type.size();j++){
							if(dist_type.get(i).get(j)>max){
								max=dist_type.get(i).get(j);
								i_max=i;
								j_max=j;
							}
						}
					}
					if(i_max==-1||j_max==-1){
						break;
					}
					
					cluster_type.get(i_max).addAll(cluster_type.get(j_max));
					cluster_type.remove(j_max);
					int[] x=new int[cluster_type_sum.get(i_max).length];
					for (int j=0;j<cluster_type_sum.get(i_max).length;j++){
						x[j]=cluster_type_sum.get(i_max)[j]+cluster_type_sum.get(j_max)[j];
					}
					cluster_type_sum.set(i_max, x);
					cluster_type_sum.remove(j_max);
					
					dist_type.remove(j_max);
					for (int i=0;i<dist_type.size();i++){
						dist_type.get(i).remove(j_max);
					}
					//ArrayList<Integer> combi_cluster=new ArrayList<Integer>();
					//combi_cluster.addAll(cluster_type.get(i_max));
					//combi_cluster.addAll(cluster_type.get(j_max));
					for (int i=0;i<dist_type.size();i++){
						if(i==i_max){
							dist_type.get(i).set(i_max,0.0);
							dist_type.get(i_max).set(i,0.0);
						}
						else{
							double d=distance_type(cluster_type_sum.get(i),cluster_type_sum.get(i_max));
							dist_type.get(i).set(i_max,d);
							dist_type.get(i_max).set(i,d);
							//dist_type.get(i).set(i_max,distance_type(cluster_type.get(i),combi_cluster));
							//dist_type.get(i_max).set(i,distance_type(cluster_type.get(i),combi_cluster));
						}
						
					}
				}
				
				
				
				ArrayList<ArrayList<Integer>> cluster_context=new ArrayList<ArrayList<Integer>>();
				ArrayList<int[][][]> cluster_context_sum=new ArrayList<int[][][]>();
				for (int i=0;i<affinity.size();i++){
					ArrayList<Integer> c=new ArrayList<Integer>();
					c.add(i);
					cluster_context.add(c);
					cluster_context_sum.add(clone(affinity.get(i)));
				}
				
				
				//hierarchial clustering in a similar manner to the previous loop for the context-dependent count vectors
				//note that the update is performed on 30CPUs in parallel to speed up the process 
				
				//System.out.println(entities[aa]);
				while(cluster_context.size()>0){
					int i_max=-1;
					int j_max=-1;
					double max=0;
					
					long t1=System.currentTimeMillis();
					for (int i=0;i<cluster_context.size();i++){
						for (int j=i+1;j<cluster_context.size();j++){
							if(dist_context.get(i).get(j)>max){
								max=dist_context.get(i).get(j);
								i_max=i;
								j_max=j;
							}
						}
					}
					long t2=System.currentTimeMillis();
					if(i_max==-1||j_max==-1){
						break;
					}
					//System.out.println("context	"+cluster_context.size()+"	"+max);
					
					cluster_context.get(i_max).addAll(cluster_context.get(j_max));
					cluster_context.remove(j_max);
					int[][][] x=new int[cluster_context_sum.get(i_max).length][cluster_context_sum.get(i_max)[0].length][cluster_context_sum.get(i_max)[0][0].length];
					for (int j1=0;j1<cluster_context_sum.get(i_max).length;j1++){
						for (int j2=0;j2<cluster_context_sum.get(i_max)[j1].length;j2++){
							for (int j3=0;j3<cluster_context_sum.get(i_max)[j1][j2].length;j3++){
								x[j1][j2][j3]=cluster_context_sum.get(i_max)[j1][j2][j3]+cluster_context_sum.get(j_max)[j1][j2][j3];
							}
						}
					}
					cluster_context_sum.set(i_max, x);
					cluster_context_sum.remove(j_max);
					
					dist_context.remove(j_max);
					for (int i=0;i<dist_context.size();i++){
						dist_context.get(i).remove(j_max);
					}
					//ArrayList<Integer> combi_cluster=new ArrayList<Integer>();
					//combi_cluster.addAll(cluster_type.get(i_max));
					//combi_cluster.addAll(cluster_type.get(j_max));
					for (int i=0;i<dist_context.size();i++){
						if(i==i_max){
							dist_context.get(i).set(i_max,0.0);
							dist_context.get(i_max).set(i,0.0);
						}
						else{
							double d=distance_context(cluster_context_sum.get(i),cluster_context_sum.get(i_max));
							dist_context.get(i).set(i_max,d);
							dist_context.get(i_max).set(i,d);
							//dist_type.get(i).set(i_max,distance_type(cluster_type.get(i),combi_cluster));
							//dist_type.get(i_max).set(i,distance_type(cluster_type.get(i),combi_cluster));
						}
					}
					
					/*
					for (int i=0;i<cluster_context.size();i++){
						vector_context[i]=distance_context(sum2(affinity,cluster_context.get(i)),sum_combi_cluster);
					}
					
					
					ArrayList<Integer> combi_cluster=new ArrayList<Integer>();
					combi_cluster.addAll(cluster_context.get(i_max));
					combi_cluster.addAll(cluster_context.get(j_max));
					long t3=System.currentTimeMillis();
					
					distance_context(combi_cluster);
					long t4=System.currentTimeMillis();
					for (int i=0;i<dist_context.size();i++){
						if(i==i_max){
							dist_context.get(i).set(i_max,0.0);
							dist_context.get(i_max).set(i,0.0);
						}
						else{
							dist_context.get(i).set(i_max,vector_context[i]);
							dist_context.get(i_max).set(i,vector_context[i]);
						}
						
					}
					dist_context.remove(j_max);
					long t5=System.currentTimeMillis();
					for (int i=0;i<dist_context.size();i++){
						dist_context.get(i).remove(j_max);
					}
					
					cluster_context.get(i_max).addAll(cluster_context.get(j_max));
					cluster_context.remove(j_max);
					long t6=System.currentTimeMillis();
					
					System.out.println(cluster_context.size()+"	"+System.currentTimeMillis()+"	"+(t2-t1)+"	"+(t3-t2)+"	"+(t4-t3)+"	"+(t5-t4)+"	"+(t6-t5));
					*/
				}
				
				
				
				int[][] count=new int[cluster_type.size()][cluster_context.size()];
				for (int i=0;i<cluster_type.size();i++){
					for (int j=0;j<cluster_type.get(i).size();j++){
						count[i][index(cluster_type.get(i).get(j),cluster_context)]++;
					}
				}
				
				
				//combine the clusters derived from type and context vectors, respectively
				//to clusters integrating both metrics
				
				//determine in which combination of type and context clusters there are samples
				ArrayList<int[]> pairs=new ArrayList<int[]>();
				for (int i=0;i<count.length;i++){
					for (int j=0;j<count[i].length;j++){
						if(count[i][j]>=1){
							pairs.add(new int[]{i,j});
						}
					}
				}
				
				//restructue the samples into the combination clusters based on their combined metrics
				ArrayList<ArrayList<Integer>> clusters=new ArrayList<ArrayList<Integer>>();
				for (int i=0;i<pairs.size();i++){
					clusters.add(new ArrayList<Integer>());
				}
				
				for (int i=0;i<cluster_type.size();i++){
					for (int j=0;j<cluster_type.get(i).size();j++){
						int index=index(new int[]{i,index(cluster_type.get(i).get(j),cluster_context)},pairs);
						if(index!=-1){
							clusters.get(index).add(cluster_type.get(i).get(j));
						}
					}
				}
				
				
				//output of the clusters for that particular entitiy
				FileWriter out=new FileWriter(file_out+entities[aa]+".txt");
				BufferedWriter output= new BufferedWriter(out);
				for (int i=0;i<names_type.size();i++){
					output.write(names_type.get(i)+"	"+entities[aa]+"	"+index(i,clusters)+"	"+index(i,cluster_type)+"	"+index(i,cluster_context));
					output.newLine();
				}
				output.close();
				
				
			}
			//System.out.println(n1+"	"+n2+"	"+n3);
			
			
		}
		catch(Exception e){
			StackTraceElement[] aa=e.getStackTrace();
			for (int i=0;i<aa.length;i++){
				System.out.println(i+"	"+aa[i].getLineNumber());
			}
			System.out.print(e);
		}
	}
	
	/*
	public static int[] type(int[] a){
		int[] b=new int[6];
		for (int i=0;i<a.length;i++){
			b[i/16]+=a[i];
		}
		return b;
	}
	
	public static int[] type(int[] a, int x){
		int[] b=new int[16];
		for (int i=0;i<b.length;i++){
			b[i]+=a[i+x*16];
		}
		return b;
	}*/
	/*
	public static int[] type(int[] a, int x, int y){
		int[] b=new int[4];
		if(y==0){
			for (int i=0;i<16;i++){
				b[i/4]+=a[i+x*16];
			}
		}
		else if(y==1){
			for (int i=0;i<16;i++){
				b[i%4]+=a[i+x*16];
			}
		}
		
		return b;
	}*/
	
	public static int[] clone(int[] x){
		int[] y=new int[x.length];
		for (int i=0;i<x.length;i++){
			y[i]=x[i];
		}
		return y;
	}
	
	public static int[][][] clone(int[][][] x){
		int[][][] y=new int[x.length][x[0].length][x[0][0].length];
		for (int i=0;i<x.length;i++){
			for (int j=0;j<x[i].length;j++){
				for (int k=0;k<x[i][j].length;k++){
					y[i][j][k]=x[i][j][k];
				}
			}
		}
		return y;
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
	
	public static double median(double[][] a){
		ArrayList<Double> list=new ArrayList<Double>();
		for (int i=0;i<a.length;i++){
			for (int j=i+1;j<a[i].length;j++){
				list.add(a[i][j]);
			}
		}
		Collections.sort(list);
		return list.get((int)(list.size()*0.5));
	}
	
	public static double logGamma(double z){
		/*if(z<gamma.length){
			return gamma[z];
		}*/
		if((int)(z*10)<gamma.length){
			return gamma[(int)(z*10)];
		}
		else{
			//System.out.println("warning "+z);
			return 0.5*Math.log(2*Math.PI)+(z-0.5)*Math.log(z)-z;
		}
	}
	
	public static int index(String s, String[] t){
		for (int i=0;i<t.length;i++){
			if(t[i].equals(s)){
				return i;
			}
		}
		return -1;
	}
	
	public static int index(int[] a, ArrayList<int[]> b){
		for (int i=0;i<b.size();i++){
			boolean equal=true;
			for (int j=0;j<a.length;j++){
				if(a[j]!=b.get(i)[j]){
					equal=false;
					break;
				}
			}
			if(equal){
				return i;
			
			}
		}
		return -1;
	}
	
	//updating the distance matrix in each step of the hierarchical clustering
	//distributing the load on n cpus
	
	/*
	public static void distance_context(ArrayList<Integer> combi_cluster){
		//int n_processor=24;
		int[][][] sum_combi_cluster=sum2(affinity,combi_cluster);
		
		vector_context=new double[cluster_context.size()];
		for (int i=0;i<cluster_context.size();i++){
			vector_context[i]=distance_context(sum2(affinity,cluster_context.get(i)),sum_combi_cluster);
		}
		
	}*/
	
	
	//updating the distance matrix in each step of the hierarchical clustering
	//distributing the load on n cpus
	/*
	public static void distance_type(ArrayList<Integer> combi_cluster){
		//int n_processor=24;
		vector_type=new double[cluster_type.size()];
		for (int i=0;i<cluster_type.size();i++){
			vector_type[i]=distance_type(cluster_type.get(i),combi_cluster);
		}
		
	
	}
	*/
	
	public static int index(int a , ArrayList<ArrayList<Integer>> c){
		for (int i=0;i<c.size();i++){
			for (int j=0;j<c.get(i).size();j++){
				if(c.get(i).get(j).intValue()==a){
					return i;
				}
			}
		}
		
		return -1;
	}
	
	public static boolean contains(String s, ArrayList<String> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).equals(s)){
				return true;
			}
		}
		return false;
	}
	
	public static boolean contains(int s, ArrayList<Integer> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i)==s){
				return true;
			}
		}
		return false;
	}
	
	
	public static double[][] initialize_type(ArrayList<int[]> type){
		
		double[][] distance_type=new double[type.size()][type.size()];
		for (int c=0;c<distance_type.length;c++){
			for (int i=0;i<distance_type[c].length;i++){
				if(i==c){
					distance_type[c][i]=0;
				}
				else if(c<i){
					distance_type[c][i]=distance_type(type.get(c),type.get(i));//,sum[c],sum[i]
					//nn++;
				}
				else{
					distance_type[c][i]=distance_type[i][c];
				}
			}
		}
		return distance_type;
	}
	
	
	public static double[][] initialize_context(ArrayList<int[][][]> affinity){
				
		double[][] distance_context=new double[affinity.size()][affinity.size()];
		for (int c=0;c<distance_context.length;c++){
			for (int i=0;i<distance_context[c].length;i++){
				if(i==c){
					distance_context[c][i]=0;
				}
				else if(c<i){
					distance_context[c][i]=distance_context(affinity.get(c),affinity.get(i));
					//System.out.println(distance_context(affinity.get(c),affinity.get(i)));
					//nn++;
				}
				else{
					distance_context[c][i]=distance_context[i][c];
				
				}
			}
		}
		return distance_context;
	
	}
	
	/*
	//the time-critical is to initalize the computation of the initial distance metrics
	//to save some time this method equally distributes this on 30 CPUs. each Subthread 
	//instance is responsible to compute 1 row in the matrix. as soon as it is done it returns
	// and the cpu is loaded with another instance
	public static void initializeX(){
		Subthread[] threads=new Subthread[distance_type.length];
		for (int i=0;i<threads.length;i++){
			threads[i]=new Subthread();
			threads[i].c=i;
		}
		ArrayList<Integer> pending=new ArrayList<Integer>();
		for (int i=0;i<distance_type.length;i++){
			pending.add(i);
		}
		ArrayList<Integer> running=new ArrayList<Integer>();
		
		while(pending.size()>0){
			for (int i=running.size()-1;i>=0;i--){
				if(threads[running.get(i)].done){
					running.remove(i);
//					System.out.println(running.get(i)+" done");
				}
			}
			
			if(running.size()<30){
				running.add(pending.get(0));
				
				threads[pending.get(0)].start();
				//System.out.println(pending.get(0)+" start");
				pending.remove(0);
			}
			
		}
		
		try{
			boolean all_done=false;
			do{
				Thread.sleep(500);
				all_done=true;
				for (int i=0;i<threads.length;i++){
					if(!threads[i].done){
						all_done=false;
						break;
					}
				}
			}while(!all_done);
		}
		catch(Exception e){
			StackTraceElement[] aa=e.getStackTrace();
			for (int i=0;i<aa.length;i++){
				System.out.println(i+"	"+aa[i].getLineNumber());
			}
			System.out.println(e);
		}
		
	}
	*/
	
	//Computes 1 row for the initial type-distance matrix
	
	/*
	private static class Subthread_Type extends Thread{
		int c_start=-1;
		int c_end=-1;
		volatile boolean done=false;
		ArrayList<Integer> combi_cluster=new ArrayList<Integer>();
		
		public void run(){
			for (int c=c_start;c<=c_end;c++){
			//	System.out.println(c+"	"+c_start+"	"+c_end+"	"+vector_context.length+"	"+cluster_context.size());
				vector_type[c]=distance_type(cluster_type.get(c),combi_cluster);
				
			}
			done=true;
		}
		
	}*/

	
	/*
	//computes 1 row for the initial context-distance matrix
	private static class Subthread_Context extends Thread{
		int c_start=-1;
		int c_end=-1;
		volatile boolean done=false;
		ArrayList<Integer> combi_cluster=new ArrayList<Integer>();
		
		public void run(){
			for (int c=c_start;c<=c_end;c++){
			//	System.out.println(c+"	"+c_start+"	"+c_end+"	"+vector_context.length+"	"+cluster_context.size());
				vector_context[c]=distance_context(cluster_context.get(c),combi_cluster);
				
			}
			done=true;
		}
	}*/
	
	/*
	//compbines 
	private static class Subthread extends Thread{
		int c=-1;
		volatile boolean done=false;
		
		public void run(){
			for (int i=0;i<distance_type[c].length;i++){
				if(i==c){
					distance_type[c][i]=0;
					distance_context[c][i]=0;
				}
				else{
					distance_type[c][i]=distance_type(type.get(c),type.get(i));
					distance_context[c][i]=distance_context(affinity.get(c),affinity.get(i));
				}
			}
			done=true;
		}
	}
	*/
	public static int sum (int[] a){
		int sum=0;
		for (int i=0;i<a.length;i++){
			sum+=a[i];
		}
		return sum;
		
	}
	
	public static double sum (double[] a){
		double sum=0;
		for (int i=0;i<a.length;i++){
			sum+=a[i];
		}
		return sum;
		
	}
	
	/*
	public static double distance_type(ArrayList<Integer> cluster1, ArrayList<Integer> cluster2){
		return distance_type(sum(type,cluster1), sum(type,cluster2));
	}
	public static double distance_context(ArrayList<Integer> cluster1, ArrayList<Integer> cluster2){
		return distance_context(sum2(affinity,cluster1), sum2(affinity,cluster2));
		//return distance_context(sum(cosmic,cluster1), sum(cosmic,cluster2));
		//TODO
	}
	*/
	
	public static int[] sum(ArrayList<int[]> type, ArrayList<Integer> cluster){
		int[] sum=new int[type.get(0).length];
		for (int i=0;i<cluster.size();i++){
			for (int j=0;j<sum.length;j++){
				sum[j]+=type.get(cluster.get(i))[j];
			}
		}
		return sum;
	}
	
	public static int[][][] sum2(ArrayList<int[][][]> affinity, ArrayList<Integer> cluster){
		int[][][] sum=new int[20][6][4];
		for (int i=0;i<cluster.size();i++){
			for (int j=0;j<sum.length;j++){
				for (int k=0;k<sum[j].length;k++){
					for (int l=0;l<sum[j][k].length;l++){
						sum[j][k][l]+=affinity.get(cluster.get(i))[j][k][l];
					}
				}
			}
		}
		return sum;
	}
	
	
	/*
	public static double distance_type(int[] vv, int[]ww, int sum_vv, int sum_ww){
		
		//double a=(double)(Math.min(sum(vv),sum(ww)))/(double)(vv.length);//(double)(sum(vv)+sum(ww))/(double)(vv.length);
		double a=(double)(Math.min(sum_vv,sum_ww))/(double)(vv.length);//(double)(sum(vv)+sum(ww))/(double)(vv.length);
		//double[] xx=new double[vv.length];
		//for (int i=0;i<vv.length;i++){
		//	xx[i]=a;
		//}
		
		return log_ratio(vv,ww,a, sum_vv, sum_ww);
	}*/
	
	public static double distance_type(int[] vv, int[]ww){
		
		double a=(double)(Math.min(sum(vv),sum(ww)))/(double)(vv.length);//(double)(sum(vv)+sum(ww))/(double)(vv.length);
		//double a=(double)(Math.min(sum_vv,sum_ww))/(double)(vv.length);//(double)(sum(vv)+sum(ww))/(double)(vv.length);
		double[] xx=new double[vv.length];
		for (int i=0;i<vv.length;i++){
			xx[i]=a;
		}
		
		return log_ratio(vv,ww,xx);
	}
	
	/*
	public static double distance_context(int[][][] v, int[][][] w){//(int[] v, int[] w){//
		double sum2=0;
		for (int i=10-5;i<=9+5;i++){
			for (int j=0;j<6;j++){
				double a=(double)(Math.min(sum(v[i][j]),sum(w[i][j])));
				//a=Math.sqrt(a);
				double[] z=new double[4];
				for (int k=0;k<4;k++){
					//System.out.println(i+","+j+","+k);
					z[k]=a*ref[i][j/3][k];
				}
				if(sum(v[i][j])>0&&sum(w[i][j])>0){
					sum2+=log_ratio(v[i][j],w[i][j],z);
				}
			
			}
		}
		return sum2;
	}*/
	
	
	public static double distance_context(int[][][] v, int[][][] w){//(int[] v, int[] w){//
		double sum2=0;
		for (int i=10-5;i<=9+5;i++){
			for (int j=0;j<6;j++){
				double a=(double)(Math.min(sum(v[i][j]),sum(w[i][j])));
				//a=Math.sqrt(a);
				double[] z=new double[4];
				for (int k=0;k<4;k++){
					//System.out.println(i+","+j+","+k);
					z[k]=a*ref[i][j/3][k];
				}
				if(sum(v[i][j])>0&&sum(w[i][j])>0){
					sum2+=log_ratio(v[i][j],w[i][j],z);
				}
			
			}
		}
		return sum2;
	}
	
//	public static double distance_context(int[][][] v, int[][][] w){//(int[] v, int[] w){//
//		double sum=0;
//		for (int i=10-5;i<=9+5;i++){
//			for (int j=0;j<6;j++){
//				double a=(double)(Math.min(sum(v[i][j]),sum(w[i][j])));
//				if(sum(v[i][j])>0&&sum(w[i][j])>0){
//					sum-=log_ratio(sum(v[i][j]),sum(w[i][j]),a);
//					//System.out.println(sum(v[i][j])+"	"+sum(w[i][j])+"	"+a+"	"+log_ratio(sum(v[i][j]),sum(w[i][j]),a));
//					for (int k=0;k<4;k++){
//						if(v[i][j][k]>0&&w[i][j][k]>0){
//							//System.out.println(log_ratio(v[i][j][k],w[i][j][k],a*ref[i][j/3][k]));
//							sum+=log_ratio(v[i][j][k],w[i][j][k],a*ref[i][j/3][k]);
//							//ddd
//						}
//					}
//				}
//				
//			}
//		}
//		return sum;
//		
//		/*
//		double sum2=0;
//		for (int i=10-5;i<=9+5;i++){
//			for (int j=0;j<6;j++){
//				double a=(double)(Math.min(sum(v[i][j]),sum(w[i][j])));
//				double[] z=new double[4];
//				for (int k=0;k<4;k++){
//					z[k]=a*ref[i][j/3][k];
//				}
//				if(sum(v[i][j])>0&&sum(w[i][j])>0){
//					sum2+=log_ratio(v[i][j],w[i][j],z);
//				}
//			
//			}
//		}
//		return sum2;*/
//	}
	
	
	/*
	public static double distance_context(int[] v, int[] w){
		double sum2=0;
		for (int i=0;i<6;i++){
			int [] x_pre=new int[4];
			int [] y_pre=new int[4];
			int [] x_post=new int[4];
			int [] y_post=new int[4];
			for (int j=0;j<16;j++){
				x_pre[j/4]+=v[i*16+j];
				y_pre[j/4]+=w[i*16+j];
				x_post[j%4]+=v[i*16+j];
				y_post[j%4]+=w[i*16+j];
			}
			if(sum(x_pre)>0&&sum(y_pre)>0){
				double a=(double)(Math.min(sum(x_pre),sum(y_pre)));//(double)(x_pre.length);
				double [] z=new double[4];
				if(i<=2){
					for (int j=0;j<z.length;j++){
						z[j]=a*ref[9][0][j];
					}
				}
				else{
					for (int j=0;j<z.length;j++){
						z[j]=a*ref[9][1][j];
					}
				}
				sum2+=log_ratio(x_pre,y_pre,z);
			}
			if(sum(x_post)>0&&sum(y_post)>0){
				double a=(double)(Math.min(sum(x_post),sum(y_post)));//(double)(x_post.length);
				double [] z=new double[4];
				if(i<=2){
					for (int j=0;j<z.length;j++){
						z[j]=a*ref[10][0][j];
					}
				}
				else{
					for (int j=0;j<z.length;j++){
						z[j]=a*ref[9][1][j];
					}
				}
				
				sum2+=log_ratio(x_post,y_post,z);
			}
			
		}
		
		
		return sum2;
		
	}	*/
	
	
	
	public static double log_ratio(int[] v, int[] w, double[] x){
		/*double sum1=0;
		double sum2=0;
		for (int i=0;i<v.length;i++){
			sum1+=v[i];
			sum2+=w[i];
		}
		double sum=0;
		for (int i=0;i<v.length;i++){
			sum-=((double)(v[i])/(double)(sum1)-(double)(w[i])/(double)(sum2))*((double)(v[i])/(double)(sum1)-(double)(w[i])/(double)(sum2));
		}
		return sum;*/
		/*double sum12=0;
		double sum11=0;
		double sum22=0;
		for (int i=0;i<v.length;i++){
			sum12+=v[i]*w[i];
			sum11+=v[i]*v[i];
			sum22+=w[i]*w[i];
		}
		return Math.log(sum12/Math.sqrt(sum11*sum22));
		*/
		double ratio=-log_ratio(sum(v),sum(w),sum(x));
		
		//System.out.println(sum(v)+"	"+sum(w)+"	"+sum(x)+"	"+log_ratio(sum(v),sum(w),sum(x)));
		for (int i=0;i<v.length;i++){
			if(v[i]>0&&w[i]>0){
				ratio+=log_ratio(v[i],w[i],x[i]);
				
			}
			//System.out.println(v[i]+"	"+w[i]+"	"+x[i]+"	"+log_ratio(v[i],w[i],x[i]));
		}
		return ratio;
	}
	
	/*
	public static double log_ratio(int[] v, int[] w, double a, int sum_v, int sum_w){
		
		double ratio=-log_ratio(sum_v,sum_w,v.length*a);//sum(x)
		
		for (int i=0;i<v.length;i++){
			ratio+=log_ratio(v[i],w[i],a);
		}
		return ratio;
	}
	*/
	
//	static int n1=0;
//	static int n2=0;
//	static int n3=0;
	public static double log_ratio(int v, int w, double x){
		/*
		if(v==0||w==0){
			n1++;
		}
		else if(v==1||w==1){
			n2++;
		}
		else{
			n3++;
		}*/
		//System.out.println(v+"	"+w+"	"+x);
		//System.out.println(logGamma(w+x)+"	"+Gamma.logGamma(w+x));
		/*
		if(v==0||w==0){
			//System.out.println(0+"	"+(logGamma(v+w+x)+logGamma(x)-logGamma(w+x)-logGamma(v+x)));
			return 0;
		}
		if(v==1){
			//System.out.println(Math.log(w+x)-Math.log(x)+"	"+(logGamma(v+w+x)+logGamma(x)-logGamma(w+x)-logGamma(v+x)));
			return Math.log(w+x)-Math.log(x);
		}
		if(v==1){
			//System.out.println(Math.log(v+x)-Math.log(x)+"	"+(logGamma(v+w+x)+logGamma(x)-logGamma(w+x)-logGamma(v+x)));
			return Math.log(v+x)-Math.log(x);
		}
		*/
		return logGamma(v+w+x)+logGamma(x)-logGamma(w+x)-logGamma(v+x);
	}
}
