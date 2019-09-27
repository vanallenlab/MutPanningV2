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
 *			separately for pan cancer. Starting from the	*
 *			clusters generated on the cancer types, this	*
 *			scripts finishes the clusters for the pan		*
 *			cancer cohort. 									*
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

import org.apache.commons.math3.special.Gamma;



public class ClusteringPanCancer {

	static String[] entities=new String[0];
	static String file_clusters="";
	static String file_type="";
	static String file_cosmic="";
	
	static String file_affinity="";
	static String file_reference="";
	
	static String file_out_cosmic="";
	static String file_out_affinity="";
	static String file_out_samples="";

	static String[] index_header_samples={"ID","Sample","Cohort"};
	
	/*
	
	
	static ArrayList<ArrayList<Integer>> cluster_type=new ArrayList<ArrayList<Integer>>();
	static ArrayList<ArrayList<Integer>> cluster_context=new ArrayList<ArrayList<Integer>>();
	static ArrayList<ArrayList<String>> entities_type=new ArrayList<ArrayList<String>>();
	static ArrayList<ArrayList<String>> entities_context=new ArrayList<ArrayList<String>>();
	
	
	static double[][] distance_type=new double[0][0];
	static double[][] distance_context=new double[0][0];
	
	
	
	static double[] vector_context=new double[0];
	*/
	
	static double[] gamma=new double[50000+1];
	static double[][][] ref=new double[20][2][4];
	
	
	/*
	 * argument0: root file
	 *  argument1: sample file (needed to extract the names of the entities only)
	 *	argument2: min no. samples per cluster
	 * argument3: min no. mutations per cluster
	 * argument4: no. processors to distribute the workload
	 */
	public static void main (String[] args){
		int threshold_samples=Integer.parseInt(args[2]);
		int threshold_mutations=Integer.parseInt(args[3]);
		
		file_clusters=args[0]+"Clustering/";
		file_type=args[0]+"AffinityCounts/TypeCount.txt";
		file_cosmic=args[0]+"AffinityCounts/CosmicCount.txt";
		
		file_affinity=args[0]+"AffinityCounts/AffinityCount.txt";
		file_reference=args[4]+"FileReferenceCount.txt";
		
		file_out_cosmic=args[0]+"ClusteringComplete/ClusteringComplete_Cosmic.txt";
		file_out_affinity=args[0]+"ClusteringComplete/ClusteringComplete_Affinity.txt";
		file_out_samples=args[0]+"ClusteringComplete/ClusteringComplete_Samples.txt";
		
		if(!new File(args[0]+"ClusteringComplete/").exists()){
			new File(args[0]+"ClusteringComplete/").mkdir();
		}
		
		//Determine all entity names, as clustering is performed separately for each Step this is needed to coordinate the order
		try{
			FileInputStream in=new FileInputStream(args[1]);
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header=index_header(input.readLine().split("	"),index_header_samples) ;
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
		
		
		//save the Gamma values as they are needed frequently later, to save some time
		for (int i=0;i<gamma.length;i++){
			gamma[i]=Gamma.logGamma((double)(i)/10.0);
		}
		
		try{
			String s="";
			
			//counts in the reference genome,needed as prior in the distance metrics
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
			
			
			
			ArrayList<String> names=new ArrayList<String>();
			ArrayList<String> clusters_type=new ArrayList<String>();
			ArrayList<String> clusters_both=new ArrayList<String>();
			ArrayList<String> clusters_context=new ArrayList<String>();
			ArrayList<String> clusters_all_type=new ArrayList<String>();
			ArrayList<String> clusters_all_both=new ArrayList<String>();
			ArrayList<String> clusters_all_context=new ArrayList<String>();
			ArrayList<String> ent=new ArrayList<String>();
			
			
			//read in the clusters from the previous clustering procedure
			for (int i=0;i<entities.length;i++){
				in=new FileInputStream(file_clusters+entities[i]+".txt");
				inn=new DataInputStream(in);
				input= new BufferedReader(new InputStreamReader(inn));
				
				while((s=input.readLine())!=null){
					names.add(s.split("	")[0]);
					ent.add(entities[i]);
					clusters_both.add(s.split("	")[1]+s.split("	")[2]);
					clusters_type.add(s.split("	")[1]+s.split("	")[3]);
					clusters_context.add(s.split("	")[1]+s.split("	")[4]);
					if(!contains(s.split("	")[1]+s.split("	")[2],clusters_all_both)){
						clusters_all_both.add(s.split("	")[1]+s.split("	")[2]);
					}
					if(!contains(s.split("	")[1]+s.split("	")[3],clusters_all_type)){
						clusters_all_type.add(s.split("	")[1]+s.split("	")[3]);
					}
					if(!contains(s.split("	")[1]+s.split("	")[4],clusters_all_context)){
						clusters_all_context.add(s.split("	")[1]+s.split("	")[4]);
					}
				}
				input.close();
			}
			
			ArrayList<String> names_type=new ArrayList<String>();
			ArrayList<int[]> type=new ArrayList<int[]>();
			ArrayList<int[]> cosmic=new ArrayList<int[]>();
			ArrayList<int[][][]> affinity=new ArrayList<int[][][]>();
			
			//read in the count vectors for the samples
			in=new FileInputStream(file_type);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			input.readLine();
			
			FileInputStream in2=new FileInputStream(file_affinity);
			DataInputStream inn2=new DataInputStream(in2);
			BufferedReader input2= new BufferedReader(new InputStreamReader(inn2));
			input2.readLine();
			
			FileInputStream in3=new FileInputStream(file_cosmic);
			DataInputStream inn3=new DataInputStream(in3);
			BufferedReader input3= new BufferedReader(new InputStreamReader(inn3));
			input3.readLine();
			
			
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				String[] t2=input2.readLine().split("	");
				String[] t3=input3.readLine().split("	");
				
				if(!t[0].equals(t2[0])){
					int a=1/0;
				}
				if(!t[0].equals(t3[0])){
					int a=1/0;
				}
				
					
					int[] c=new int[6];
					for (int i=0;i<c.length;i++){
						c[i]=Integer.parseInt(t[i+1]);
					}
					int[] c3=new int[96];
					for (int i=0;i<c3.length;i++){
						c3[i]=Integer.parseInt(t3[i+1]);
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
						cosmic.add(c3);
						affinity.add(c2);
					}
					
				
			}
			input.close();
			input2.close();
			
			
			ArrayList<ArrayList<Integer>> cluster_type=new ArrayList<ArrayList<Integer>>();
			ArrayList<ArrayList<Integer>> cluster_context=new ArrayList<ArrayList<Integer>>();
			ArrayList<int[]> cluster_type_sum=new ArrayList<int[]>();
			ArrayList<int[][][]> cluster_context_sum=new ArrayList<int[][][]>();
			
			ArrayList<ArrayList<String>> entities_type=new ArrayList<ArrayList<String>>();
			ArrayList<ArrayList<String>> entities_context=new ArrayList<ArrayList<String>>();
			
			//for (int i=0;i<clusters_all_type.size();i++){
			for (int i=0;i<clusters_all_both.size();i++){
				cluster_type.add(new ArrayList<Integer>());
				entities_type.add(new ArrayList<String>());
			}
			//for (int i=0;i<clusters_all_context.size();i++){
			for (int i=0;i<clusters_all_both.size();i++){
				cluster_context.add(new ArrayList<Integer>());
				entities_context.add(new ArrayList<String>());
			}
			
			
			// merge the entities only if they dont overlap!!
			for (int i=0;i<names_type.size();i++){
				cluster_type.get(index(clusters_both.get(index(names_type.get(i),names)),clusters_all_both)).add(i);
				cluster_context.get(index(clusters_both.get(index(names_type.get(i),names)),clusters_all_both)).add(i);
				
				if(!contains(ent.get(index(names_type.get(i),names)),entities_type.get(index(clusters_both.get(index(names_type.get(i),names)),clusters_all_both)))){
					entities_type.get(index(clusters_both.get(index(names_type.get(i),names)),clusters_all_both)).add(ent.get(index(names_type.get(i),names)));
				}
				if(!contains(ent.get(index(names_type.get(i),names)),entities_context.get(index(clusters_both.get(index(names_type.get(i),names)),clusters_all_both)))){
					entities_context.get(index(clusters_both.get(index(names_type.get(i),names)),clusters_all_both)).add(ent.get(index(names_type.get(i),names)));
				}
				
			}
			
			for (int i=0;i<cluster_type.size();i++){
				int[] x=new int[6];
				for (int j=0;j<cluster_type.get(i).size();j++){
					for (int k=0;k<type.get(cluster_type.get(i).get(j)).length;k++){
						x[k]+=type.get(cluster_type.get(i).get(j))[k];
					}
				}
				cluster_type_sum.add(x);
			}
			
			for (int i=0;i<cluster_context.size();i++){
				int[][][] x=new int[20][6][4];
				for (int j=0;j<cluster_context.get(i).size();j++){
					for (int k1=0;k1<affinity.get(cluster_context.get(i).get(j)).length;k1++){
						for (int k2=0;k2<affinity.get(cluster_context.get(i).get(j))[k1].length;k2++){
							for (int k3=0;k3<affinity.get(cluster_context.get(i).get(j))[k1][k2].length;k3++){
								x[k1][k2][k3]+=affinity.get(cluster_context.get(i).get(j))[k1][k2][k3];
							}
						}
					}
				}
				cluster_context_sum.add(x);
			}
			
			System.out.println("XXXX");
			//double[][] distance_type=new double[cluster_type.size()][cluster_type.size()];
			//double[][] distance_context=new double[cluster_context.size()][cluster_context.size()];
			//in parallel to the previous script initialization of the distance matrix in parallel
			double[][] distance_type=initialize_type(cluster_type_sum);
			System.out.println("XXXX");
			double[][] distance_context=initialize_context(cluster_context_sum);
			System.out.println("XXXX");
			
			
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
			
			//hierarchical clustering based on mutation type counts
			while(cluster_type.size()>0){
				int i_max=-1;
				int j_max=-1;
				double max=0;//median_type;//-100000000;
				for (int i=0;i<cluster_type.size();i++){
					for (int j=i+1;j<cluster_type.size();j++){
						if(overlap(entities_type.get(i),entities_type.get(j)).size()==0&&dist_type.get(i).get(j)>max){
							max=dist_type.get(i).get(j);
							i_max=i;
							j_max=j;
						}
					}
				}
				//System.out.println("type	"+cluster_type.size()+"	"+max);
				if(i_max==-1||j_max==-1){
					break;
				}
				
				//ArrayList<Integer> combi_cluster=new ArrayList<Integer>();
				//combi_cluster.addAll(cluster_type.get(i_max));
				//combi_cluster.addAll(cluster_type.get(j_max));
				
				dist_type.remove(j_max);
				for (int i=0;i<dist_type.size();i++){
					dist_type.get(i).remove(j_max);
				}
				cluster_type.get(i_max).addAll(cluster_type.get(j_max));
				cluster_type.remove(j_max);
				
				int[] x=new int[cluster_type_sum.get(i_max).length];
				for (int j=0;j<cluster_type_sum.get(i_max).length;j++){
					x[j]=cluster_type_sum.get(i_max)[j]+cluster_type_sum.get(j_max)[j];
				}
				cluster_type_sum.set(i_max, x);
				cluster_type_sum.remove(j_max);
				
				for (int i=0;i<dist_type.size();i++){
					if(i==i_max){
						dist_type.get(i).set(i_max,0.0);
						dist_type.get(i_max).set(i,0.0);
					}
					else{
						double d=distance_type(cluster_type_sum.get(i),cluster_type_sum.get(i_max));
						dist_type.get(i).set(i_max,d);
						dist_type.get(i_max).set(i,d);
					}
				}
				
				entities_type.get(i_max).addAll(entities_type.get(j_max));
				entities_type.remove(j_max);
			}
			
			
			
			//hierarchical clustering based on mutation context counts
			while(cluster_context.size()>0){
				int i_max=-1;
				int j_max=-1;
				double max=0;
				for (int i=0;i<cluster_context.size();i++){
					for (int j=i+1;j<cluster_context.size();j++){
						if(overlap(entities_context.get(i),entities_context.get(j)).size()==0&&dist_context.get(i).get(j)>max){
							max=dist_context.get(i).get(j);
							i_max=i;
							j_max=j;
						}
					}
				}
				//System.out.println("context	"+cluster_context.size()+"	"+max);
				if(i_max==-1||j_max==-1){
					break;
				}
				
				//ArrayList<Integer> combi_cluster=new ArrayList<Integer>();
				//combi_cluster.addAll(cluster_context.get(i_max));
				//combi_cluster.addAll(cluster_context.get(j_max));
				
				//distance_context(combi_cluster);
				
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
				*/
				
				entities_context.get(i_max).addAll(entities_context.get(j_max));
				entities_context.remove(j_max);
			}
			
			
			//extract all the pairs of clusters that contain at least 1 sample
			int[][] count=new int[cluster_type.size()][cluster_context.size()];
			for (int i=0;i<cluster_type.size();i++){
				for (int j=0;j<cluster_type.get(i).size();j++){
					count[i][index(cluster_type.get(i).get(j),cluster_context)]++;
				}
			}
			
			ArrayList<int[]> pairs=new ArrayList<int[]>();
			for (int i=0;i<count.length;i++){
				for (int j=0;j<count[i].length;j++){
					if(count[i][j]>=1){
						pairs.add(new int[]{i,j});
					}
				}
			}
			
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
			
			//small clusters are not considered (they might lead to overfitting and blow up the
			//resources for the subsequent steps. standard criteria to throw a cluster:
			// - less than 3 samples
			// - less than 1000 mutations
			
			for (int i=clusters.size()-1;i>=0;i--){
				if(clusters.get(i).size()<threshold_samples||sum(sum(type,clusters.get(i)))<threshold_mutations){
					clusters.remove(i);
				}
			}
			
			ArrayList<int[]> cluster_ty=new ArrayList<int[]>();
			
			
			FileWriter out=new FileWriter(file_out_cosmic);
			BufferedWriter output= new BufferedWriter(out);
			output.write("ClusterID");
			for (int i=1;i<=96;i++){
				output.write("	Bin_"+i);
			}
			output.newLine();
			
			for (int i=0;i<clusters.size();i++){
				int[] sum=sum(cosmic,clusters.get(i));
				output.write(i+"");
				for (int j=0;j<sum.length;j++){
					output.write("	"+sum[j]);
				}
				output.newLine();
				
				int[] sum2=sum(type,clusters.get(i));
				cluster_ty.add(sum2);
			}
			output.close();
			
			out=new FileWriter(file_out_affinity);
			output= new BufferedWriter(out);
			output.write("ClusterID");
			for (int i=0;i<20;i++){
				for (int j=0;j<6;j++){
					for (int k=0;k<4;k++){
						int p=0;
						if(i<10){
							p=i-10;
						}
						else{
							p=i-9;
						}
						output.write("	Position_"+p+"_Type_"+(j+1)+"_"+new String[]{"A","C","G","T"}[k]);
					}
				}
			}
			output.newLine();
			
			ArrayList<int[][][]> cluster_affinity=new ArrayList<int[][][]>();
			for (int i=0;i<clusters.size();i++){
				int[][][] sum=sum2(affinity,clusters.get(i));
				
				output.write(i+"");
				for (int j=0;j<sum.length;j++){
					for (int k=0;k<sum[j].length;k++){
						for (int l=0;l<sum[j][k].length;l++){
							output.write("	"+sum[j][k][l]);
						}
					}
				}
				output.newLine();
				cluster_affinity.add(sum);
			}
			output.close();
			
			out=new FileWriter(file_out_samples);
			output= new BufferedWriter(out);
			output.write("SampleID	Entity	ClusterID");
			output.newLine();
			for (int i=0;i<affinity.size();i++){
				int c=index(i,clusters);
				if(c==-1){
					c=max_likelihood(affinity.get(i), type.get(i), cluster_affinity, cluster_ty);
				}
				output.write(names_type.get(i)+"	"+ent.get(index(names_type.get(i),names))+"	"+c);
				output.newLine();
			}
			output.close();
			
			
			
		}
		catch(Exception e){
			StackTraceElement[] aa=e.getStackTrace();
			for (int i=0;i<aa.length;i++){
				System.out.println(i+"	"+aa[i].getLineNumber());
			}
			
			System.out.println(e);
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
	
	//for samples that do not fall into any cluster determine the cluster  that they would best fit in
	public static int max_likelihood(int[][][] affinity, int[] cosmic, ArrayList<int[][][]> cluster_affinity, ArrayList<int[]> cluster_cosmic){
		double max=-10000000;
		int i_max=-1;
		for (int i=0;i<cluster_affinity.size();i++){
			double prob=distance_type(cosmic,cluster_cosmic.get(i))+distance_context(affinity,cluster_affinity.get(i));
			if(prob>max){
				max=prob;
				i_max=i;
			}
		}
		return i_max;
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
	
	//method to make sure that not clusters within the same entity are merged together
	public static ArrayList<String> overlap(ArrayList<String> list1,ArrayList<String> list2){
		ArrayList<String> overlap=new ArrayList<String>();
		for (int i=0;i<list1.size();i++){
			for (int j=0;j<list2.size();j++){
				if(list1.get(i).equals(list2.get(j))){
					if(!contains(list2.get(j),overlap)){
						overlap.add(list2.get(j));
					}
				}
			}
		}
		return overlap;
	}
	
	//equally distribute the computational load to update the distance matrix on 30 CPUs
	
	/*
	static void distance_context(ArrayList<Integer> combi_cluster){
		//int n_processor=30;
		vector_context=new double[cluster_context.size()];
		
		int a=cluster_context.size()/n_processor+1;
		Subthread_Context[] threads=new Subthread_Context[n_processor];
		for (int i=0;i<threads.length;i++){
			threads[i]=new Subthread_Context();
			threads[i].combi_cluster=combi_cluster;
			threads[i].c_start=i*a;
			threads[i].c_end=Math.min(cluster_context.size()-1, i*a+a-1);
			threads[i].start();
		}
		
		
		try{
			boolean all_done=false;
			do{
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
			System.out.println(e);
		}
	}*/
	
	/*
	private static class Subthread_Context extends Thread{
		int c_start=-1;
		int c_end=-1;
		volatile boolean done=false;
		ArrayList<Integer> combi_cluster=new ArrayList<Integer>();
		
		public void run(){
			for (int c=c_start;c<=c_end;c++){
				vector_context[c]=distance_context(cluster_context.get(c),combi_cluster);
			}
			done=true;
		}
	}*/
	
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
	
	public static int index(String s, ArrayList<String> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).equals(s)){
				return i;
			}
		}
		return -1;
	}
	
/*
	public static void initialize(){
		
		Subthread[] threads=new Subthread[cluster_type.size()];
		for (int i=0;i<threads.length;i++){
			threads[i]=new Subthread();
			threads[i].c=i;
		}
		ArrayList<Integer> pending=new ArrayList<Integer>();
		for (int i=0;i<cluster_type.size();i++){
			pending.add(i);
		}
		ArrayList<Integer> running=new ArrayList<Integer>();
		
		while(pending.size()>0){
			for (int i=running.size()-1;i>=0;i--){
				if(threads[running.get(i)].done){
					running.remove(i);
				}
			}
			
			if(running.size()<30){
				running.add(pending.get(0));
				
				threads[pending.get(0)].start();
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
		
	}*/
	
	public static boolean contains(String s, ArrayList<String> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).equals(s)){
				return true;
			}
		}
		return false;
	}
	
	//thread responsible to calculate 1 row of the distance matrices, used in parallelization step 
	
	
	public static double[][] initialize_type(ArrayList<int[]> cluster_type){
		double[][] distance_type=new double[cluster_type.size()][cluster_type.size()];
		for (int c=0;c<cluster_type.size();c++){
			for (int i=0;i<distance_type[c].length;i++){
				if(i==c){
					distance_type[c][i]=0;
				}
				else if(c<i){
					distance_type[c][i]=distance_type(cluster_type.get(c),cluster_type.get(i));
				}
				else{
					distance_type[c][i]=distance_type[i][c];
				}
			}
		}
		return distance_type;
	}
	
	public static double[][] initialize_context(ArrayList<int[][][]> cluster_context){
		double[][] distance_context=new double[cluster_context.size()][cluster_context.size()];
		for (int c=0;c<cluster_context.size();c++){
			for (int i=0;i<distance_context[c].length;i++){
				if(i==c){
					distance_context[c][i]=0;
				}
				else if(c<i){
					distance_context[c][i]=distance_context(cluster_context.get(c),cluster_context.get(i));
				}
				else{
					distance_context[c][i]=distance_context[i][c];
				}
			}
		}
		return distance_context;
	}
	
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
	public static double distance_type(ArrayList<Integer> cluster1, ArrayList<Integer> cluster2){
		return distance_type(sum(type,cluster1), sum(type,cluster2));
	}
	public static double distance_context(ArrayList<Integer> cluster1, ArrayList<Integer> cluster2){
		return distance_context(sum2(affinity,cluster1), sum2(affinity,cluster2));
	}*/
	
	public static double distance_type(int[] vv, int[] ww){
		
		double a=Math.sqrt((double)(Math.min(sum(vv),sum(ww)))/(double)(vv.length));//);//);//(double)(sum(vv)+sum(ww))/(double)(vv.length);
		double[] xx=new double[vv.length];
		for (int i=0;i<vv.length;i++){
			xx[i]=a;
		}
		double sum1=log_ratio(vv,ww,xx);
		return sum1;
	}
	
	public static double distance_typeX(int[] vv, int[] ww){
		
		double a=(double)(Math.min(sum(vv),sum(ww)))/(double)(vv.length);//(double)(sum(vv)+sum(ww))/(double)(vv.length);
		double[] xx=new double[vv.length];
		for (int i=0;i<vv.length;i++){
			xx[i]=a;
		}
		double sum1=log_ratio(vv,ww,xx);
		return sum1;
	}
	
	
	public static double distance_contextX(int[][][] v, int[][][] w){//(int[] v, int[] w){//
		double sum2=0;
		for (int i=10-5;i<=9+5;i++){
			for (int j=0;j<6;j++){
				double a=(double)(Math.min(sum(v[i][j]),sum(w[i][j])));
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
	
	
	public static double distance_context(int[][][] v, int[][][] w){//(int[] v, int[] w){//
		double sum2=0;
		for (int i=10-5;i<=9+5;i++){
			for (int j=0;j<6;j++){
				double a=Math.sqrt((double)(Math.min(sum(v[i][j]),sum(w[i][j]))));//);//);
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
	
	public static double log_ratio(int[] v, int[] w, double[] x){
		double ratio=-log_ratio(sum(v),sum(w),sum(x));
		
		for (int i=0;i<v.length;i++){
			ratio+=log_ratio(v[i],w[i],x[i]);
		}
		return ratio;
	}
	public static double log_ratio(int v, int w, double x){
		return logGamma(v+w+x)+logGamma(x)-logGamma(w+x)-logGamma(v+x);
	}

	public static double logGamma(double z){
		
		if((int)(z*10)<gamma.length){
			return gamma[(int)(z*10)];
		}
		else{
			return 0.5*Math.log(2*Math.PI)+(z-0.5)*Math.log(z)-z;
		}
	}
}
