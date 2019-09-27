/************************************************************           
 * MutPanning - Main Class for Desktop Verion				*
 * 															*   
 * Author:		Felix Dietlein								*   
 *															*   
 * Copyright:	(C) 2019 									*   
 *															*   
 * License:		BSD-3-Clause open source license			*   
 *															*   
 * Summary: MutPanning is detects cancer genes, based on	*
 * their mutation counts as well as the sequence context 	*
 * around mutations. This script coordinates the execution	*
 * of all substeps. All arguments are facultative except	*
 * of the root file in which the script should write its 	*
 * (intermediate) output files. If no other arguments are 	*
 * provided the script expects the maf file, the sample 	*
 * file and the Hg19 folder in this folder. Otherwise, paths*
 * to maf file, sample file and Hg19 folder can be provided *
 * as separate arguments. Please launch this script with    *
 * -Xmx10G in the VM arguments to allocate enough memory. 	*
 * Designed to run on a Windows, Linux or MacOS platform with*
 * 1 CPU, 8 GB RAM and 3GB hard disk space (2.75GB for 	*
 * Hg19 folder and 350 MB for intermediate output files.	*
 *************************************************************/


import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Cursor;
import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
//import java.io.FileFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Collections;

import javax.imageio.ImageIO;
import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JRadioButton;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.UIManager;
import javax.swing.border.Border;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;

public class MutPanningDesktop {
	//static MutPanningDesktop x=new MutPanningDesktop();
	//public static int ii=0;
	//static Frame frame=null;
	//static JButton button_continue=null;
	//static JButton button_back=null;
	
	
	
	static int page_intro=0;
	static int page_sample=1;
	static int page_maf=2;
	static int page_hg19=3;
	static int page_intermediate=4;
	static int page_parameters=5;
	static int page_review=6;
	static int page_run=7;
	static int page_completion=8;
	
	
	static String[] index_header_maf={"Hugo_Symbol","Chromosome","Start_Position","End_Position","Strand","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","Tumor_Sample_Barcode"};
	static String[] index_header_samples_short={"Sample","Cohort"};
	static String[] index_header_samples={"ID","Sample","Cohort"};
	
	static String[] file_checklist={
		"Exons_Hg19.txt",
		"FileReferenceCount.txt",
		"TableBessel.txt",
		"SignificanceFilter/11.ooc",
		"SignificanceFilter/ApprovedSymbols.txt",
		"SignificanceFilter/BlackList.txt",
		"SignificanceFilter/BlackListOLD.txt",
		"SignificanceFilter/blat_mac",
		"SignificanceFilter/blat_linux",
		"SignificanceFilter/census.csv",
		"SignificanceFilter/hg19.2bit",
		"AnnotationHg19/Annotation_chr1.txt",
		"AnnotationHg19/Annotation_chr2.txt",
		"AnnotationHg19/Annotation_chr3.txt",
		"AnnotationHg19/Annotation_chr4.txt",
		"AnnotationHg19/Annotation_chr5.txt",
		"AnnotationHg19/Annotation_chr6.txt",
		"AnnotationHg19/Annotation_chr7.txt",
		"AnnotationHg19/Annotation_chr8.txt",
		"AnnotationHg19/Annotation_chr9.txt",
		"AnnotationHg19/Annotation_chr10.txt",
		"AnnotationHg19/Annotation_chr11.txt",
		"AnnotationHg19/Annotation_chr12.txt",
		"AnnotationHg19/Annotation_chr13.txt",
		"AnnotationHg19/Annotation_chr14.txt",
		"AnnotationHg19/Annotation_chr15.txt",
		"AnnotationHg19/Annotation_chr16.txt",
		"AnnotationHg19/Annotation_chr17.txt",
		"AnnotationHg19/Annotation_chr18.txt",
		"AnnotationHg19/Annotation_chr19.txt",
		"AnnotationHg19/Annotation_chr20.txt",
		"AnnotationHg19/Annotation_chr21.txt",
		"AnnotationHg19/Annotation_chr22.txt",
		"AnnotationHg19/Annotation_chrX.txt",
		"AnnotationHg19/Annotation_chrY.txt",
	};
	
	public static void main(String[] args){
//		File[] ff=new File("./").listFiles();
//		for (int i=0;i<ff.length;i++){
//			System.out.println(ff[i].getAbsolutePath());
//		}
//		System.out.println(new File("./BannerLogosXV2.png").exists());
//		System.exit(0);
		
		/*
		try {
			
			
		    //System.setOut(new PrintStream(new File("Q:\\output-file.txt")));
		} catch (Exception e) {
		     e.printStackTrace();
		}
		System.out.println("TS");
		System.out.println("TS");
		System.out.println("TS");
		
		System.exit(0);
		*/
		try{
			Frame frame=new Frame(0);
			ClassLoader classLoader = Thread.currentThread().getContextClassLoader();
			frame.setIconImage(new ImageIcon(classLoader.getResource("MutPanningSymbol_Icon.png")).getImage());
			//frame.setDefaultCloseOperation();
			//frame.setSize(800,500-25);
			//frame.setVisible(true);
		}
		catch(Exception e){
			System.out.println(e);
		}
	}
	
	private static class Frame extends JFrame implements ActionListener, KeyListener, MouseListener{
		boolean is_running=false;
		boolean success=false;
		
		JLabel label_top=null;
		JLabel[] label_text=null;
		JLabel label_link1=null; 
		JLabel label_link2=null;
		public  int ii=0;
		JButton button_continue=null;
		JButton button_back=null;
		Font f=new Font("Helvetica",Font.PLAIN,12);
		Font f2=new Font("Helvetica",Font.BOLD,14);
		JPanel panel=null;
		JPanel pp=null;
		JLabel[] label_side=null;
		Bullet[] bullets_side=null;
		
		JButton select_file=null;
		JTextField input_path=null;
		
		
		JRadioButton option1=null;
		JRadioButton option2=null;
		ButtonGroup group=null;
		JLabel[] label_param=null;
		JTextField[] input_param=null;
		
		String path_sample_file="";//"S:\\TestRunsMutPanning\\SamplesSkin.txt";//"S:\\TestRunMutPanning2\\SamplesComplete.txt";
		String path_maf_file="";//"S:\\TestRunsMutPanning\\MutationsSkin.maf";//"S:\\TestRunMutPanning2\\MutationsComplete.maf";
		String path_hg19_file="";//"Q:\\Hg19";//"S:\\TestRunMutPanning2\\Hg19\\";
		String path_temporary_file="";//"S:\\TestRunMutPanning3";
		String[] params={"3","1000","100","5000"};
		boolean standard=true;
		
		String status="";
		int progress=0;
		JProgressBar progress_bar=null;
		
		Subthread thread=new Subthread();
		public void reset(){
			//panel.removeAll();
			pp.removeAll();
			if(label_top!=null){
				label_top.setVisible(false);
				//label_top=null;
			}
			if(label_text!=null){
				for (int i=0;i<label_text.length;i++){
					label_text[i].setVisible(false);
					label_text[i]=null;
				}
				label_text=null;
			}
			if(label_link1!=null){
				label_link1.setVisible(false);
				label_link1=null;
			}
			if(label_link2!=null){
				label_link2.setVisible(false);
				label_link2=null;
			}
			button_back.setEnabled(true);
			button_continue.setEnabled(true);
			
			if(select_file!=null){
				select_file.setVisible(false);
				select_file=null;
			}
			
			if(input_path!=null){
				input_path.setVisible(false);
				input_path=null;
			}
			
			if(option1!=null){
				option1.setVisible(false);
				option1=null;
			}
			if(option2!=null){
				option2.setVisible(false);
				option2=null;
			}
			if(group!=null){
				group=null;
			}
			if(label_param!=null){
				for (int i=0;i<label_param.length;i++){
					label_param[i].setVisible(false);
					label_param[i]=null;
				}
				label_param=null;
			}
			if(input_param!=null){
				for (int i=0;i<input_param.length;i++){
					input_param[i].setVisible(false);
					input_param[i]=null;
				}
				input_param=null;
			}
			
			
		}
		
		public Frame(int ii){//, MutPanningDesktop ref
			this.ii=ii;
			//System.out.println(ii);
			this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			this.setSize(800,500-25);
			this.setTitle("MutPanning");
			//super();
			setLayout(new FlowLayout());
			 panel = (JPanel) getContentPane();
			
			panel.setLayout(null);
			
			label_side=new JLabel[]{new JLabel("Introduction"),new JLabel("Select Sample File"), new JLabel("Select Maf File"), new JLabel("Path to Hg19 Folder"), new JLabel("Path to Temporary Files"),new JLabel("Set Parameters"), new JLabel("Review"), new JLabel("Run MutPanning"), new JLabel("Completion")};
			for (int i=0;i<label_side.length;i++){
				label_side[i].setFont(f2);
				label_side[i].setBounds(40, 50+30*i, label_side[i].preferredSize().width,  label_side[i].preferredSize().height);
				label_side[i].setForeground(Color.GRAY);
				label_side[i].addMouseListener(this);
				panel.add(label_side[i]);
			}
			bullets_side=new Bullet[label_side.length];
			String OS = System.getProperty("os.name").toLowerCase();
			for (int i=0;i<bullets_side.length;i++){
				if(i==ii){
					bullets_side[i]=new Bullet(true);
				}
				else{
					bullets_side[i]=new Bullet(false);
				}
				if(OS.indexOf("win") >= 0){
					bullets_side[i].setBounds(15, 55+30*i, 10, 10);	
				}
				else {
					bullets_side[i].setBounds(15, 55+30*i-4, 10, 10);	
				}
				bullets_side[i].addMouseListener(this);
				panel.add(bullets_side[i]);
			}
			
			ClassLoader classLoader = Thread.currentThread().getContextClassLoader();
			//classLoader.getResource(BannerLogosXV2.png)
			
			//InputStream input = classLoader.getResourceAsStream("test.png");
			//new JLabel(ImageIO.read(input));
			//JLabel ll=new JLabel(new ImageIcon(new ImageIcon("resources/MutPanningSymbolSmallXV2.png").getImage().getScaledInstance(150,91, Image.SCALE_SMOOTH)));
			JLabel ll=new JLabel(new ImageIcon(new ImageIcon(classLoader.getResource("MutPanningSymbolSmallXV2.png")).getImage().getScaledInstance(150,91, Image.SCALE_SMOOTH)));
			ll.setBounds(-10, 270, 200, 200);
			panel.add(ll);
			
			//JLabel lll=new JLabel(new ImageIcon(new ImageIcon("Q:\\BannerLogosXV2.png").getImage().getScaledInstance(496,45, Image.SCALE_SMOOTH)));
			JLabel lll=new JLabel(new ImageIcon(new ImageIcon(classLoader.getResource("BannerLogosXV2.png")).getImage().getScaledInstance(330,30, Image.SCALE_SMOOTH)));
			lll.setBounds(40, 270+30, 900, 200);
			panel.add(lll);
			
			button_back=new JButton("Go Back");
			button_back.addActionListener(this);
			button_continue=new JButton("Continue");
			button_continue.addActionListener(this);
			button_continue.setBounds(220+555-150,330+15,150,30);
			button_back.setBounds(220+555-350,330+15,150,30);
			
			panel.add(button_back);
			panel.add(button_continue);
			
			pp=new JPanel();
			pp.setLayout(null);
			pp.setBounds(220, 50, 555, 280);
			pp.setBackground(Color.WHITE);
			pp.setBorder(BorderFactory.createLineBorder(Color.black));//new Border(.createEmptyBorder()
			panel.add(pp);
			
			paint(); 
			
			this.setVisible(true);
		}
		
		
		public boolean review_sample(){
			File f=new File(path_sample_file);
			if(!f.exists()){
				JOptionPane.showMessageDialog(this, "File "+ path_sample_file +" does not exist.");
				return false;
			}
			if(!f.isFile()){
				JOptionPane.showMessageDialog(this, "File "+ path_sample_file +" is not a file.");
				return false;
			}
			if(!f.canRead()){
				JOptionPane.showMessageDialog(this, "Insufficient permissions. Cannot read file "+path_sample_file);
				return false;
			}
			
			try{
				FileInputStream in=new FileInputStream(path_sample_file);
				DataInputStream inn=new DataInputStream(in);
				BufferedReader input= new BufferedReader(new InputStreamReader(inn));
				String[] t=input.readLine().split("	");
				input.close();
				for (int i=0;i<index_header_samples_short.length;i++){
					if(index(index_header_samples_short[i],t)==-1){
						JOptionPane.showMessageDialog(this, "Column header do not match. Column \""+index_header_samples_short[i]+"\" not found.");
						return false;
					}
				}
			}
			catch(Exception e){ 
				JOptionPane.showMessageDialog(this, "Could not open "+path_sample_file);
				return false;
			}
			
			return true;
		}
		public boolean review_maf(){
			File f=new File(path_maf_file);
			if(!f.exists()){
				JOptionPane.showMessageDialog(this, "File "+ path_hg19_file +" does not exist.");
				return false;
			}
			if(!f.isFile()){
				JOptionPane.showMessageDialog(this, "File "+ path_hg19_file +" is not a file.");
				return false;
			}
			if(!f.canRead()){
				JOptionPane.showMessageDialog(this, "Insufficient permissions. Cannot read file + "+path_maf_file);
				return false;
			}
			try{
				FileInputStream in=new FileInputStream(path_maf_file);
				DataInputStream inn=new DataInputStream(in);
				BufferedReader input= new BufferedReader(new InputStreamReader(inn));
				String[] t=input.readLine().split("	");
				input.close();
				for (int i=0;i<index_header_maf.length;i++){
					if(index(index_header_maf[i],t)==-1){
						JOptionPane.showMessageDialog(this, "Column header do not match. Column \""+index_header_maf[i]+"\" not found in "+path_maf_file);
						return false;
					}
				}
				
			}
			catch(Exception e){
				JOptionPane.showMessageDialog(this, "Could not open "+path_sample_file);
				return false;
			}
			
			return true;
		}
		
		public static int index(String s, String[] t){
			for (int i=0;i<t.length;i++){
				if(t[i].equals(s)){
					return i;
				}
			}
			return -1;
		}
		
		public boolean review_hg19(){
			if(!path_hg19_file.substring(path_hg19_file.length()-1, path_hg19_file.length()).equals("/")&&!path_hg19_file.substring(path_hg19_file.length()-1, path_hg19_file.length()).equals("\\")){
				path_hg19_file=path_hg19_file+"/";
			}
			File f=new File(path_hg19_file);
			if(!f.exists()){
				JOptionPane.showMessageDialog(this, "File "+ path_hg19_file +" does not exist.");
				return false;
			}
			if(!f.isDirectory()){
				JOptionPane.showMessageDialog(this, "File "+ path_hg19_file +" is not a directory.");
				return false;
			}
			if(!f.canRead()){
				JOptionPane.showMessageDialog(this, "Insufficient permissions. Cannot read folder "+ path_hg19_file);
				return false;
			}
			
			
			
			for (int i=0;i<file_checklist.length;i++){
				if(!new File(path_hg19_file+file_checklist[i]).exists()){
					JOptionPane.showMessageDialog(this, "Hg19 folder is not complete. "+file_checklist[i]+" not found. Please download the latest version.");
					return false;
				}
			}
			return true;
		}
		
		public boolean review_intermediate(){
			
			if(!path_temporary_file.substring(path_temporary_file.length()-1, path_temporary_file.length()).equals("/")&&!path_temporary_file.substring(path_temporary_file.length()-1, path_temporary_file.length()).equals("\\")){
				path_temporary_file=path_temporary_file+"/";
			}
			
			File f=new File(path_temporary_file);
			if(!f.exists()){
				JOptionPane.showMessageDialog(this, "File "+ path_temporary_file +" does not exist.");
				return false;
			}
			if(!f.isDirectory()){
				JOptionPane.showMessageDialog(this, "File "+ path_temporary_file +" is not a directory.");
				return false;
			}
			if(!f.canRead()){
				JOptionPane.showMessageDialog(this, "Insufficient permissions. Cannot read folder "+ path_temporary_file);
				return false;
			}
			if(!f.canWrite()){
				JOptionPane.showMessageDialog(this, "Insufficient permissions. Cannot write to folder "+ path_temporary_file);
				return false;
			}
			
			return true;
		}
		public boolean review_parameters(){
			if(standard){
				return true;
			}
			else{
				for (int i=0;i<params.length;i++){
					if(!isInteger(params[i])){
						JOptionPane.showMessageDialog(this, "All parameters need to be integer numbers.");
						return false;
					}
					else{
						if(Integer.parseInt(params[i])<0){
							JOptionPane.showMessageDialog(this, "All parameters need to be non-negative.");
							return false;
						}
					}
				}
				return true;
			}
		}

		public static boolean isInteger(String s){
			try{
				Integer.parseInt(s);
				return true;
			}
			catch(Exception e){
				return false;
			}
		}
		
		public boolean sanity_check(){
			
			if(ii==page_sample){
				return review_sample();
			}
			else if(ii==page_maf){
				return review_maf();
			}
			else if(ii==page_hg19){
				return review_hg19();
			}
			else if(ii==page_intermediate){
				return review_intermediate();
			}
			else if(ii==page_parameters){
				return review_parameters();
			}
			else if(ii==page_review){
				return review_sample()&&review_maf()&&review_hg19()&&review_intermediate()&&review_parameters();
			}
			return true;	
		}
		
		public void paint (){
			for (int i=0;i<label_side.length;i++){
				label_side[i].setForeground(Color.GRAY);
				bullets_side[i].active=false;
			}
			
			label_side[ii].setForeground(Color.BLUE);
			bullets_side[ii].active=true;
			if(ii==page_intro){
			
				button_back.setEnabled(false);
				button_continue.setText("Continue");
				
				label_top=new JLabel("Welcome to MutPanning");
				label_top.setFont(f);
				label_top.setBounds(220+10, 30, label_top.preferredSize().width, label_top.preferredSize().height);
				label_top.setForeground(Color.BLACK);
				label_top.setVisible(true);
				panel.add(label_top);
				
				String[] text={
						"MutPanning analyzes aggregated DNA sequencing data of tumor patients to identify genes that are",
						"likely to be functionally relevant, based on their abundance of nonsynonymous mutations or their",
						"increased number of mutations in unusual nucleotide contexts that deviate from the background ",
						"mutational process. You can run the algorithm for multiple cancer types at the same time.",
						"To run MutPanning, you need the following three files:",
						"- <b>Mutation File</b> (*.maf): mutations that you would like to include in the analysis",
						"- <b>Sample File</b> (*.txt): contains sample IDs and assocates them with cancer types",
						"- <b>Hg19 Folder</b>: contains the Hg19 annotations, such as nucleotide contexts (1.1 GB disk space).",
						"If you have not downloaded it already, you can download this folder",
						"- <b>empty folder</b>: this folder will be used to write temporary files and MutPanning results",
						"<b>System requirements</b>: 8GB memory, 2 GHz CPU, 2 cores, 400MB free disk space.",
						"This guide will let you select these files on your computer and then launch MutPanning locally.",
						"If your computer does not meet these criteria, you can run MutPanning online on "//<a href=\"genepattern.org\">genepattern.org</a>
				};
				int[] offset={
					0,
					0,
					0,
					0,
					0,
					20,
					20,
					20,
					100,
					20,
					0,
					0,
					0
				};
				
				//Windows, Linux or MacOSX system,
				label_text=new JLabel[text.length];
				for (int i=0;i<text.length;i++){
					label_text[i]=new JLabel("<html>"+text[i]+"</html>");
					label_text[i].setFont(f);
					label_text[i].setBounds(10+offset[i],10+20*i,label_text[i].getPreferredSize().width,label_text[i].getPreferredSize().height);
					label_text[i].setVisible(true);
					pp.add(label_text[i]);
				}
				
				label_link1=new JLabel("<html>"+"<a href=\"http://storage.googleapis.com/mutpanning_hg19/Hg19.zip\">here</a>"+"</html>");
				label_link1.setFont(f);
				label_link1.setBounds(478,10+20*8,label_link1.getPreferredSize().width,label_link1.getPreferredSize().height);
				goWebsite(label_link1,"http://storage.googleapis.com/mutpanning_hg19/Hg19.zip");
				label_link1.setVisible(true);
				pp.add(label_link1);
				label_link1.setCursor(new Cursor(Cursor.HAND_CURSOR));
				
				label_link2=new JLabel("<html>"+"<a href=\"genepattern.org\">genepattern.org</a>"+"</html>");
				label_link2.setFont(f);
				label_link2.setBounds(448,10+20*12,label_link2.getPreferredSize().width,label_link2.getPreferredSize().height);
				goWebsite(label_link2,"http://software.broadinstitute.org/cancer/software/genepattern/quick-start");
				label_link2.setVisible(true);
				pp.add(label_link2);
				label_link2.setCursor(new Cursor(Cursor.HAND_CURSOR));
				
			}
			else if(ii==page_sample){
				button_continue.setText("Continue");
				label_top=new JLabel("MutPanning - Select Sample File");
				label_top.setFont(f);
				label_top.setBounds(220+10, 30, label_top.preferredSize().width, label_top.preferredSize().height);
				label_top.setForeground(Color.BLACK);
				label_top.setVisible(true);
				panel.add(label_top);
				String[] text={
					"Please select a <b>sample file</b>. The sample file should be a tab-delimited *.txt file that contains",
					"at least two columns labeled <b>Sample</b> and <b>Cohort</b>. The Sample column contains unique",
					"sample name for each sample in the maf file. No duplicates allowed, all samples in the mutation ",
					"file must be listed in the sample file. The cohort column contains the subcohort (e.g. cancer type) ",
					"to which the sample belongs (case-sensitive).",
					"If you are unsure about the file format, you can download an exemplary sample file "
				};
				
				label_text=new JLabel[text.length];
				for (int i=0;i<text.length;i++){
					label_text[i]=new JLabel("<html>"+text[i]+"</html>");
					label_text[i].setFont(f);
					label_text[i].setBounds(10,10+20*i,label_text[i].getPreferredSize().width,label_text[i].getPreferredSize().height);
					label_text[i].setVisible(true);
					pp.add(label_text[i]);
				}
				
				label_link1=new JLabel("<html>"+"<a href=\"http://storage.googleapis.com/mutpanning_hg19/SamplesSkin.txt\">here</a>"+"</html>");
				label_link1.setFont(f);
				label_link1.setBounds(460,10+20*5,label_link1.getPreferredSize().width,label_link1.getPreferredSize().height);
				goWebsite(label_link1,"http://storage.googleapis.com/mutpanning_hg19/SamplesSkin.txt");
				label_link1.setVisible(true);
				pp.add(label_link1);
				label_link1.setCursor(new Cursor(Cursor.HAND_CURSOR));
				
				label_link2=new JLabel("<html>"+"<b>Path to sample file</b>"+"</html>");
				label_link2.setFont(f);
				label_link2.setBounds(30,10+20*7,label_link2.getPreferredSize().width,label_link2.getPreferredSize().height);
				label_link2.setVisible(true);
				pp.add(label_link2);
				
				input_path=new JTextField();
				input_path.setFont(f);
				input_path.setText(path_sample_file);
				input_path.setBounds(30,10+20*8,300,25);
				input_path.setVisible(true);
				input_path.addActionListener(this);
				input_path.addKeyListener(this);
				pp.add(input_path);
				select_file=new JButton("...");
				select_file.setBounds(350,10+20*8,50,25);
				select_file.addActionListener(this);
				select_file.setVisible(true);
				pp.add(select_file);
				

			
			}
			else if(ii==page_maf){
				button_continue.setText("Continue");
				label_top=new JLabel("MutPanning - Select Maf File");
				label_top.setFont(f);
				label_top.setBounds(220+10, 30, label_top.preferredSize().width, label_top.preferredSize().height);
				label_top.setForeground(Color.BLACK);
				label_top.setVisible(true);
				panel.add(label_top);
				String[] text={
					"Please select a <b>mutation file</b>. The mutation file should be a tab-delimited stanard *.maf format",
					"and contain the following columns: <b>Hugo_Symbol</b> (gene name), <b>Chromosome</b>, <b>Start_Position</b>",
					"(according to Hg19), <b>End_Position</b> (according to Hg19), <b>Strand</b> (1 or -1), <b>Variant_Classification</b>",
					" (e.g. Missense_Mutation, Nonsense_Mutation, In_Frame_Del, Silent, etc.), <b>Variant_Type</b> (e.g. ",
					"SNP, DEL, INS, etc.), <b>Reference_Allele</b> (reference nucleotide in Hg19), <b>Tumor_Seq_Allele1</b> ",
					"(allele A in tumor), <b>Tumor_Seq_Allele2</b> (allele B in tumor), <b>Tumor_Sample_Barcode</b>.",
					"Names in column <b>Tumor_Sample_Barcode</b> should exactly match the names in the sample file ",
					"(case-sensitive). If you are unsure about the file format, you can download an example ",
					"Make sure that <b>Variant_Classification</b> and <b>Variant_Type</b> are annotated by exactly the same",
					"terms/names as in the exemplary mutation file (case-senitive, hyphenation, underscores, etc.)"
				};
				
				label_text=new JLabel[text.length];
				for (int i=0;i<text.length;i++){
					label_text[i]=new JLabel("<html>"+text[i]+"</html>");
					label_text[i].setFont(f);
					label_text[i].setBounds(10,10+20*i,label_text[i].getPreferredSize().width,label_text[i].getPreferredSize().height);
					label_text[i].setVisible(true);
					pp.add(label_text[i]);
				}
				
				label_link1=new JLabel("<html>"+"<a href=\"http://storage.googleapis.com/mutpanning_hg19/MutationsSkin.maf\">here</a>"+"</html>");
				label_link1.setFont(f);
				label_link1.setBounds(482,10+20*7,label_link1.getPreferredSize().width,label_link1.getPreferredSize().height);
				goWebsite(label_link1,"http://storage.googleapis.com/mutpanning_hg19/MutationsSkin.maf");
				label_link1.setVisible(true);
				pp.add(label_link1);
				label_link1.setCursor(new Cursor(Cursor.HAND_CURSOR));
				
				
				label_link2=new JLabel("<html>"+"<b>Path to mutation file</b>"+"</html>");
				label_link2.setFont(f);
				label_link2.setBounds(30,10+20*10+10,label_link2.getPreferredSize().width,label_link2.getPreferredSize().height);
				label_link2.setVisible(true);
				pp.add(label_link2);
				
				
				input_path=new JTextField();
				input_path.setFont(f);
				input_path.setText(path_maf_file);
				input_path.setBounds(30,10+20*11+10,300,25);
				input_path.setVisible(true);
				input_path.addActionListener(this);
				input_path.addKeyListener(this);
				pp.add(input_path);
				select_file=new JButton("...");
				select_file.setBounds(350,10+20*11+10,50,25);
				select_file.addActionListener(this);
				select_file.setVisible(true);
				pp.add(select_file);
				
			}
			else if(ii==page_hg19){
				button_continue.setText("Continue");
				label_top=new JLabel("MutPanning - Select Hg19 Folder");
				label_top.setFont(f);
				label_top.setBounds(220+10, 30, label_top.preferredSize().width, label_top.preferredSize().height);
				label_top.setForeground(Color.BLACK);
				label_top.setVisible(true);
				panel.add(label_top);
				String[] text={
					"Please select the path to the <b>Hg19 folder</b>. This folder contains all the information that ",
					"MutPanning needs about the reference genome, e.g. to determine the nucleotide contexts around ",
					"mutations. If you have not downloaded this folder already, you can download it ",
					"1.1GB free disk space required. Please make sure to unzip Hg19.zip before you proceed. ",
					
				};
				
				label_text=new JLabel[text.length];
				for (int i=0;i<text.length;i++){
					label_text[i]=new JLabel("<html>"+text[i]+"</html>");
					label_text[i].setFont(f);
					label_text[i].setBounds(10,10+20*i,label_text[i].getPreferredSize().width,label_text[i].getPreferredSize().height);
					label_text[i].setVisible(true);
					pp.add(label_text[i]);
				}
				
				
				label_link1=new JLabel("<html>"+"<a href=\"http://storage.googleapis.com/mutpanning_hg19/Hg19.zip\">here</a>"+"</html>");
				label_link1.setFont(f);
				label_link1.setBounds(438,10+20*2,label_link1.getPreferredSize().width,label_link1.getPreferredSize().height);
				goWebsite(label_link1,"http://storage.googleapis.com/mutpanning_hg19/Hg19.zip");
				label_link1.setVisible(true);
				pp.add(label_link1);
				label_link1.setCursor(new Cursor(Cursor.HAND_CURSOR));
				
				
				label_link2=new JLabel("<html>"+"<b>Path to Hg19 folder</b>"+"</html>");
				label_link2.setFont(f);
				label_link2.setBounds(30,10+20*5,label_link2.getPreferredSize().width,label_link2.getPreferredSize().height);
				label_link2.setVisible(true);
				pp.add(label_link2);
				
				
				input_path=new JTextField();
				input_path.setFont(f);
				input_path.setText(path_hg19_file);
				input_path.setBounds(30,10+20*6,300,25);
				input_path.setVisible(true);
				input_path.addActionListener(this);
				input_path.addKeyListener(this);
				pp.add(input_path);
				select_file=new JButton("...");
				select_file.setBounds(350,10+20*6,50,25);
				select_file.addActionListener(this);
				select_file.setVisible(true);
				pp.add(select_file);
				
			}
			else if(ii==page_intermediate){
				button_continue.setText("Continue");
				label_top=new JLabel("MutPanning - Select Folder to Temporary Files");
				label_top.setFont(f);
				label_top.setBounds(220+10, 30, label_top.preferredSize().width, label_top.preferredSize().height);
				label_top.setForeground(Color.BLACK);
				label_top.setVisible(true);
				panel.add(label_top);
				String[] text={
					"Please select the path to a folder for <b>temporary files</b>. This folder will be used to write",
					"temporary files while running MutPanning. This folder will also contain the final result files ",
					"after completion of MutPanning. Ideally, this folder should be <b>completely empty</b> when launching",
					"MutPanning. If you use a folder that was used for previous MutPanning runs, MutPanning",
					"may overwrite results from previous runs. Further, please make sure that the disk has at least",
					"<b>400 MB free space</b>, to have enough space for intermediate files while running MutPanning.",
					
				};
				
				label_text=new JLabel[text.length];
				for (int i=0;i<text.length;i++){
					label_text[i]=new JLabel("<html>"+text[i]+"</html>");
					label_text[i].setFont(f);
					label_text[i].setBounds(10,10+20*i,label_text[i].getPreferredSize().width,label_text[i].getPreferredSize().height);
					label_text[i].setVisible(true);
					pp.add(label_text[i]);
				}
				
				label_link2=new JLabel("<html>"+"<b>Path to Temporary Files</b>"+"</html>");
				label_link2.setFont(f);
				label_link2.setBounds(30,10+20*7,label_link2.getPreferredSize().width,label_link2.getPreferredSize().height);
				label_link2.setVisible(true);
				pp.add(label_link2);
				
				
				input_path=new JTextField();
				input_path.setFont(f);
				input_path.setText(path_temporary_file);
				input_path.setBounds(30,10+20*8,300,25);
				input_path.setVisible(true);
				input_path.addActionListener(this);
				input_path.addKeyListener(this);
				pp.add(input_path);
				select_file=new JButton("...");
				select_file.setBounds(350,10+20*8,50,25);
				select_file.addActionListener(this);
				select_file.setVisible(true);
				pp.add(select_file);
				
			}
			else if(ii==page_parameters){
				button_continue.setText("Continue");
				label_top=new JLabel("MutPanning - Select Parameters");
				label_top.setFont(f);
				label_top.setBounds(220+10, 30, label_top.preferredSize().width, label_top.preferredSize().height);
				label_top.setForeground(Color.BLACK);
				label_top.setVisible(true);
				panel.add(label_top);
				String[] text={
					"In this step you can manually change additional <b>parameters</b> to launch MutPanning. ",
					"Unless you are familiar with the MutPanning algorithm, we recommend running MutPanning ",
					"with standard parameters to ensure stability and completion.",
					
				};
				
				label_text=new JLabel[text.length];
				for (int i=0;i<text.length;i++){
					label_text[i]=new JLabel("<html>"+text[i]+"</html>");
					label_text[i].setFont(f);
					label_text[i].setBounds(10,10+20*i,label_text[i].getPreferredSize().width,label_text[i].getPreferredSize().height);
					label_text[i].setVisible(true);
					pp.add(label_text[i]);
				}
				
				
				
				option1=new JRadioButton("<html>Run with standard parameters <b>(recommended)</b></html>");
				option1.setBackground(Color.WHITE);
				option1.setFont(f);
				option1.setBounds(30,10+20*4+10,option1.getPreferredSize().width,option1.getPreferredSize().height);
				option1.setVisible(true);
				option2=new JRadioButton("<html>Run with individual parameters <b>(not recommended)</b></html>");
				option2.setBackground(Color.WHITE);
				option2.setFont(f);
				option2.setBounds(30,10+20*6,option2.getPreferredSize().width,option2.getPreferredSize().height);
				option2.setVisible(true);
				pp.add(option1);
				pp.add(option2);
				group = new ButtonGroup();
				if(standard){
					option1.setSelected(true);
					option2.setSelected(false);
				}
				else{
					option2.setSelected(true);
					option1.setSelected(false);
				}
				group.add(option1);
				group.add(option2);
				option1.addActionListener(this);
				option2.addActionListener(this);
				String[] text_param={"minimal no. samples per cluster","minimal no. mutations per cluster","min no. samples for Bayesian model","min no. mutations for Bayesian model"};
				
				label_param=new JLabel[text_param.length];
				input_param=new JTextField[text_param.length];
				for (int i=0;i<text_param.length;i++){
					label_param[i]=new JLabel("<html><b>"+text_param[i]+"</b></html>");
					label_param[i].setFont(f);
					label_param[i].setBounds(50+10,10+7*20+10+25*i+2,label_param[i].getPreferredSize().width,label_param[i].getPreferredSize().height);
					label_param[i].setVisible(true);
					label_param[i].setEnabled(!standard);
					pp.add(label_param[i]);
					
					input_param[i]=new JTextField(""+params[i]);
					input_param[i].setFont(f);
					input_param[i].setBounds(270+10,10+7*20+10+25*i,50,label_param[i].getPreferredSize().height+5);
					input_param[i].setVisible(true);
					input_param[i].setEnabled(!standard);
					input_param[i].addKeyListener(this);
					pp.add(input_param[i]);
				}
			}
			else if(ii==page_review){
				button_continue.setText("Run MutPanning");
				label_top=new JLabel("MutPanning - Review");
				label_top.setFont(f);
				label_top.setBounds(220+10, 30, label_top.preferredSize().width, label_top.preferredSize().height);
				label_top.setForeground(Color.BLACK);
				label_top.setVisible(true);
				panel.add(label_top);
				String[] text={
					"You are about to run MutPanning with the parameters below. Please review the parameters ",
					"and click on <b>Run MutPanning</b> when you are ready. Go back, if you would like to make any changes.",
					"",
					"<b>Sample File: </b>",
					"<b>Mutation File: </b>",
					"<b>Hg19 Folder: </b>",
					"<b>Temporary Files: </b>",
					"<b>Parameters: </b>",
					path_sample_file,
					path_maf_file,
					path_hg19_file,
					path_temporary_file,
					""
				};
				if(standard){
					text[text.length-1]="standard configuration";
				}
				else{
					text[text.length-1]="non-standard configuration ("+params[0]+", "+params[1]+", "+params[2]+", "+params[3]+")";
				}
				
				
				int[] offset={
					0,
					0,
					0,
					30,
					30,
					30,
					30,
					30
				};
				
				label_text=new JLabel[text.length];
				for (int i=0;i<text.length;i++){
					if(i<offset.length){
						label_text[i]=new JLabel("<html>"+text[i]+"</html>");
						label_text[i].setFont(f);
						label_text[i].setBounds(10+offset[i],10+20*i,label_text[i].getPreferredSize().width,label_text[i].getPreferredSize().height);
						label_text[i].setVisible(true);
						pp.add(label_text[i]);
					}
					else{
						label_text[i]=new JLabel("<html>"+text[i]+"</html>");
						label_text[i].setFont(f);
						label_text[i].setBounds(150,10+20*(i-offset.length+3),label_text[i].getPreferredSize().width,label_text[i].getPreferredSize().height);
						label_text[i].setVisible(true);
						pp.add(label_text[i]);
					
					}
				}
				
				
			}
			else if(ii==page_run){
				//System.out.println("XXXXXX");
				button_continue.setText("Continue");
				button_continue.setEnabled(false);
				button_back.setEnabled(false);
				label_top=new JLabel("MutPanning - Running");
				label_top.setFont(f);
				label_top.setBounds(220+10, 30, label_top.preferredSize().width, label_top.preferredSize().height);
				label_top.setForeground(Color.BLACK);
				label_top.setVisible(true);
				panel.add(label_top);
				String[] text={
					"<b>MutPanning is running.</b>",
					"Depending on the job size and the local computer, a MutPanning run typically takes 20 minutes ",
					"to 3 hours. You can follow the progress below. Please do not close this window while MutPanning ",
					"is running, since this would cause the algorithm to stop.",
					
				};
				label_text=new JLabel[text.length];
				for (int i=0;i<text.length;i++){
					label_text[i]=new JLabel("<html>"+text[i]+"</html>");
					label_text[i].setFont(f);
					label_text[i].setBounds(10,10+20*i,label_text[i].getPreferredSize().width,label_text[i].getPreferredSize().height);
					label_text[i].setVisible(true);
					pp.add(label_text[i]);
				}
				
				if(!status.equals("")){
					label_link2=new JLabel("<html>"+"<b>Step "+progress+"/14: "+status+"</b>"+"</html>");
					label_link2.setFont(f2);
					label_link2.setBounds(30,10+20*5,label_link2.getPreferredSize().width,label_link2.getPreferredSize().height);
					label_link2.setVisible(true);
					//if(progress_bar==null){
						progress_bar=new JProgressBar(SwingConstants.HORIZONTAL, 0,14);
					//}
						progress_bar.setBounds(30, 10+20*6+10, 400, 25);
						progress_bar.setVisible(true);
						progress_bar.setValue(progress);
					pp.add(progress_bar);
					pp.add(label_link2);
				}
				
			}
			else if(ii==page_completion){
				String output_path="";
				if(new File(path_temporary_file).exists()&&new File(path_temporary_file).isDirectory()){
					File[] ff=new File(path_temporary_file).listFiles();
					for (int i=0;i<ff.length;i++){
						//System.out.println(ff[i].getAbsolutePath());
						if(ff[i].isDirectory()&&ff[i].getAbsolutePath().substring(path_temporary_file.length()).contains("SignificanceFiltered")){
							output_path=ff[i].getAbsolutePath();
						}
					}
				}
				//System.out.println(path_temporary_file);
				//System.out.println(success+"	"+output_path);
				String[] text={};
				if(success&&!output_path.equals("")){
					text=new String[]{
							"<b>MutPanning completed successfully.</b>",
							"Results are available at",
						};
				}
				else{
					text=new String[]{
							"<b>MutPanning terminated with errors</b>",
							"We apologize for the inconvenience. If this was a large job or a computer with small memory, ",
							"we recommend that you try the online version on GenePattern or the Dockerized version that ",
							"can be run in the cloud or cluster. Memory limitations can cause MutPanning to run ouf of memory.",
							"A log file is available at "//+path_temporary_file+"Log.txt"
						};
				}
				
				button_continue.setText("Close");
				button_back.setVisible(false);
				label_top=new JLabel("MutPanning - Completed");
				label_top.setFont(f);
				label_top.setBounds(220+10, 30, label_top.preferredSize().width, label_top.preferredSize().height);
				label_top.setForeground(Color.BLACK);
				label_top.setVisible(true);
				panel.add(label_top);
				
				label_text=new JLabel[text.length];
				for (int i=0;i<text.length;i++){
					label_text[i]=new JLabel("<html>"+text[i]+"</html>");
					label_text[i].setFont(f);
					label_text[i].setBounds(10,10+20*i,label_text[i].getPreferredSize().width,label_text[i].getPreferredSize().height);
					label_text[i].setVisible(true);
					pp.add(label_text[i]);
				}
				
				if(success&&!output_path.equals("")){
					label_link1=new JLabel("<html>"+"<a href=\""+output_path+"\">"+output_path+"</a>"+"</html>");
					label_link1.setFont(f);
					label_link1.setBounds(30,10+20*2,label_link1.getPreferredSize().width,label_link1.getPreferredSize().height);
					goFolder(label_link1,output_path);
					label_link1.setVisible(true);
					pp.add(label_link1);
					label_link1.setCursor(new Cursor(Cursor.HAND_CURSOR));
				}
				else{
					label_link1=new JLabel("<html>"+"<a href=\""+path_temporary_file+"Log.txt"+"\">"+path_temporary_file+"Log.txt"+"</a>"+"</html>");
					label_link1.setFont(f);
					label_link1.setBounds(138,10+20*4,label_link1.getPreferredSize().width,label_link1.getPreferredSize().height);
					goFolder(label_link1,path_temporary_file+"Log.txt");
					label_link1.setVisible(true);
					pp.add(label_link1);
					label_link1.setCursor(new Cursor(Cursor.HAND_CURSOR));
					
					
				}
				
				
			}
			else{
				label_top=new JLabel("MutPanningStep");
				label_top.setFont(f);
				label_top.setBounds(220+10, 30, label_top.preferredSize().width, label_top.preferredSize().height);
				label_top.setForeground(Color.BLACK);
				label_top.setVisible(true);
				panel.add(label_top);
			}
			this.repaint();
			//System.out.println("YYY");
		}
		
		@Override
		public void actionPerformed(ActionEvent arg0) {
			if(arg0.getSource()==button_back&&ii>0){
				//if(sanity_check()){
					ii--;
					this.reset();
					this.paint();
				//}
				//frame.setVisible(false);
				//frame=new Frame(ii);
			}
			else if(arg0.getSource()==button_continue&&ii==page_review){
				if(sanity_check()){
					if(!is_running){
						//System.out.println("A");
						ii=page_run;
						is_running=true;
						this.reset();
						this.paint();
						//System.out.println("A");
						launch_mutpanning();
						//success=launch_mutpanning();
						//System.out.println("A");
						//is_running=false;
						//ii=page_completion;
						//this.reset();
						//this.paint();
						//System.out.println("A");
					}
				}
				
			}
			else if(arg0.getSource()==button_continue&&ii<page_completion){
				if(sanity_check()){
					ii++;
					this.reset();
					this.paint();
				}
				
				//frame.setVisible(false);
				//frame=new Frame(ii);
			}
			else if(arg0.getSource()==button_continue&&ii==page_completion){
				System.exit(0);
			}
			else if(arg0.getSource()==select_file&&ii==page_sample){
				JFileChooser fc = new JFileChooser();
				fc.setFileFilter(new FileNameExtensionFilter("Text files", "txt"));
				int returnVal = fc.showOpenDialog(this);
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					path_sample_file=fc.getSelectedFile().getAbsolutePath();
		            input_path.setText(fc.getSelectedFile().getAbsolutePath());
		        } else {
		           
		        }
			}
			else if(arg0.getSource()==input_path&&ii==page_sample){
				path_sample_file=input_path.getText();
			}
			else if(arg0.getSource()==select_file&&ii==page_maf){
				JFileChooser fc = new JFileChooser();
				fc.setFileFilter(new FileNameExtensionFilter("Maf files", "maf"));
				int returnVal = fc.showOpenDialog(this);
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					path_maf_file=fc.getSelectedFile().getAbsolutePath();
		            input_path.setText(fc.getSelectedFile().getAbsolutePath());
		        } else {
		           
		        }
			}
			else if(arg0.getSource()==input_path&&ii==page_maf){
				path_maf_file=input_path.getText();
			}
			else if(arg0.getSource()==select_file&&ii==page_hg19){
				JFileChooser fc = new JFileChooser();
				fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				//fc.setFileFilter(new FileNameExtensionFilter("Maf files", "maf"));
				int returnVal = fc.showOpenDialog(this);
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					path_hg19_file=fc.getSelectedFile().getAbsolutePath();
		            input_path.setText(fc.getSelectedFile().getAbsolutePath());
		        } else {
		           
		        }
			}
			else if(arg0.getSource()==input_path&&ii==page_hg19){
				path_hg19_file=input_path.getText();
			}
			else if(arg0.getSource()==select_file&&ii==page_intermediate){
				
				/*
				JFrame frame = new JFrame();

	            JFileChooser fc = new JFileChooser();

	            fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
	            fc.setFileFilter( new FileFilter(){

	                @Override
	                public boolean accept(File f) {
	                    return f.isDirectory();
	                }

	                @Override
	                public String getDescription() {
	                    return "Any folder";
	                }

	            });

	            fc.setDialogType(JFileChooser.SAVE_DIALOG);
	            fc.setApproveButtonText("Select");

	            frame.getContentPane().add(fc);
	            frame.setVisible(true);
				
	            ArrayList<JPanel> jpanels = new ArrayList<JPanel>();

	            for(Component c : fc.getComponents()){
	                if( c instanceof JPanel ){
	                    jpanels.add((JPanel)c);
	                }
	            }

	            jpanels.get(0).getComponent(0).setVisible(false);

	            frame.pack();
	            */
	            /*
				JFrame frame = new JFrame();
				JFileChooser fc = new JFileChooser();
				fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				FileFilter filter=new FileFilter(){
	            	@Override
	                public boolean accept(File f) {
	                    return f.isDirectory();
	                }

	                public String getDescription() {
	                    return "Any folder";
	                }
	            };
				
				fc.setFileFilter(filter);
				fc.setDialogType(JFileChooser.SAVE_DIALOG);
	            fc.setApproveButtonText("Select");
	            frame.getContentPane().add(fc);
	            frame.setVisible(true);
	            */
	            
	            
				
				JFileChooser fc = new JFileChooser();
				//disableNewFolderButton(fc);
				fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				//fc.setFileFilter(new FileNameExtensionFilter("Maf files", "maf"));
				int returnVal = fc.showOpenDialog(this);
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					path_temporary_file=fc.getSelectedFile().getAbsolutePath();
		            input_path.setText(fc.getSelectedFile().getAbsolutePath());
		        } else {
		           
		        }
			}
			else if(arg0.getSource()==input_path&&ii==page_intermediate){
				path_temporary_file=input_path.getText();
			}
			else if((arg0.getSource()==option1||arg0.getSource()==option2)&&ii==page_parameters){
				standard=option1.isSelected();
				for (int i=0;i<label_param.length;i++){
					label_param[i].setEnabled(!standard);
					input_param[i].setEnabled(!standard);
				}
				repaint();
			}
			
			//System.out.println(ii);
		}
/*
		public static void disableNewFolderButton(Container c) {
		    int len = c.getComponentCount();
		    for (int i = 0; i < len; i++) {
		      Component comp = c.getComponent(i);
		      if (comp instanceof JButton) {
		        JButton b = (JButton) comp;
		        Icon icon = b.getIcon();
		        if (icon != null
		            && icon == UIManager.getIcon("FileChooser.newFolderIcon"))
		          b.setEnabled(true);
		      } else if (comp instanceof Container) {
		        disableNewFolderButton((Container) comp);
		      }
		    }
		  }*/
		
		
		@Override
		public void keyPressed(KeyEvent arg0) {
		
		}

		@Override
		public void keyReleased(KeyEvent arg0) {
			if(arg0.getSource()==input_path&&ii==page_sample){
				path_sample_file=input_path.getText();
			}
			else if(arg0.getSource()==input_path&&ii==page_maf){
				path_maf_file=input_path.getText();
			}
			else if(arg0.getSource()==input_path&&ii==page_hg19){
				path_hg19_file=input_path.getText();
			}
			else if(arg0.getSource()==input_path&&ii==page_intermediate){
				path_temporary_file=input_path.getText();
			}
			else if(ii==page_parameters){
				if(arg0.getSource()==input_param[0]){
					params[0]=input_param[0].getText();
				}
				if(arg0.getSource()==input_param[1]){
					params[1]=input_param[1].getText();
				}
				if(arg0.getSource()==input_param[2]){
					params[2]=input_param[2].getText();
				}
				if(arg0.getSource()==input_param[3]){
					params[3]=input_param[3].getText();
				}
			}
		}

		@Override
		public void keyTyped(KeyEvent arg0) {
			
		}

		@Override
		public void mouseClicked(MouseEvent arg0) {
			
			if(ii!=page_run&&ii!=page_completion){
				for (int i=0;i<bullets_side.length;i++){
					if(arg0.getSource()==bullets_side[i]||arg0.getSource()==label_side[i]){
						if(i==page_run||i==page_completion){
							JOptionPane.showMessageDialog(this, "You cannot jump to this page directly. Other steps need to be completed previously.");
						}
						else{
							if(i>ii){
								if(sanity_check()){
									ii=i;
									this.reset();
									this.paint();
								}
							}
							else{
								ii=i;
								this.reset();
								this.paint();
							}
							
							
						}
						
					}
				}
			}
			
		
			//System.out.println("Click");
		}

		@Override
		public void mouseEntered(MouseEvent arg0) {
			
		}

		@Override
		public void mouseExited(MouseEvent arg0) {
			
		}

		@Override
		public void mousePressed(MouseEvent arg0) {
			
		}

		@Override
		public void mouseReleased(MouseEvent arg0) {
			
		}

		
		/*
		@Override
		public void actionPerformed(ActionEvent e) {
			if(e.getSource()==button_back&&ii>0){
				ii--;
				updatee();
			}
			else if(e.getSource()==button_continue){
				ii++;
				updatee();
			}
			
		}*/
		
		
		public void update_status(String status, int progress){
			if(ii==page_run){
				this.status=status;
				this.progress=progress;
				reset();
				paint();
			}
		}
		
		private class Subthread extends Thread{
			
			
			/*
			
	        try{
	        	

		        // Set file print stream.
		         if(standard){
	        		update_status("Step 1",0.1);
	 	        	Thread.sleep(2000);
	 	        	update_status("Step 1",0.1);
		 	        Thread.sleep(2000);
		 	        update_status("Step 1",0.1);
			 	    Thread.sleep(2000);	
	 	        }
	 	        else{
	 	        	
	 	        }
	        }
	        catch(Exception e){
	        	return false;
	        }
	       
	        //launch mutpanning here
	        // Set console print stream.
	       */
			
			public void run(){
				
				is_running=true;
				PrintStream ps_console = System.out;
				try{
					File file = new File(path_temporary_file+"Log.txt");
			        FileOutputStream fos = new FileOutputStream(file);
			        PrintStream ps = new PrintStream(fos);
			        System.setOut(ps);
		        	
					update_status("Launching MutPanning",0);
					String[] arguments=null;//new String[11];
					if(standard){
						arguments=new String[]{path_temporary_file,path_maf_file,path_sample_file,path_hg19_file};
					}
					else{
						arguments=new String[]{path_temporary_file,path_maf_file,path_sample_file,path_hg19_file,params[0],params[1],params[2],params[3]};
						
					}
					
					String[] args=new String[11];
					for (int i=0;i<args.length;i++){
						args[i]="";
					}
					for (int i=0;i<arguments.length;i++){
						if(!arguments[i].equals("")&&!arguments[i].equals(" ")){
							args[i]=arguments[i];
						}
					}
					
					
					if(args[0].charAt(args[0].length()-1)!=java.io.File.separatorChar){
						args[0]=args[0]+java.io.File.separatorChar;
					}
					if(args[1].equals("")){
						args[1]=args[0]+"MutationsComplete.maf";
					}
					if(args[2].equals("")){
						args[2]=args[0]+"SamplesComplete.txt";
					}
					if(args[3].equals("")){
						args[3]=args[0]+"Hg19/";
					}
					if(args[4].equals("")){
						args[4]="3";
					}
					if(args[5].equals("")){
						args[5]="1000";
					}
					if(args[6].equals("")){
						args[6]="100";
					}
					if(args[7].equals("")){
						args[7]="5000";
					}
					
				
					System.out.println("Launching MutPanning with root file "+arguments[0]);
					//System.out.println(System.currentTimeMillis());
					
					
					reindex_samples(args[2],args[0]+"SamplesCompleteReindex.txt");
					args[2]=args[0]+"SamplesCompleteReindex.txt";
					
					update_status("Aligning Maf to Hg19",1);
					System.out.println("Aligning Maf to Hg19");// - "+chr[i]);
					//System.out.println(System.currentTimeMillis());
					AlignHG19.main(new String[]{args[0],args[1],args[2],args[3]});//,chr[i]});
					
					update_status("Determining Count Vectors",2);
					System.out.println("Determining Count Vectors");
					//System.out.println(System.currentTimeMillis());
					AffinityCount.main(new String[]{args[0],args[2],args[3]});
					AffinityCount_Cosmic.main(new String[]{args[0],args[2],args[3]});
					
					update_status("Clustering for Each Cancer Type",3);
					System.out.println("Clustering for Each Cancer Type");
					//System.out.println(System.currentTimeMillis());
					ClusteringEntity.main(new String[]{args[0],args[2],args[3]});
					
					update_status("Clustering PanCancer",4);
					System.out.println("Clustering PanCancer");
					//System.out.println(System.currentTimeMillis());
					ClusteringPanCancer.main(new String[]{args[0],args[2],args[4],args[5],args[3]});

					
					
					
					update_status("Count Destructive Mutations",5);
					System.out.println("Count Destructive Mutations");
					//System.out.println(System.currentTimeMillis());
					CountDestructiveMutations.main(new String[]{args[0],args[1],args[2],args[3]});
					
					
					update_status("Compute Mutation Rate Entities/Clusters",6);
					System.out.println("Compute Mutation Rate Entities/Clusters");
					//System.out.println(System.currentTimeMillis());
					ComputeMutationRateClusters_Entities.main(new String[]{args[0],args[2],args[3]});
						
					
					
					update_status("Reformat files for Bayesian model",7);
					System.out.println("Reformat files for CBASE");
					//System.out.println(System.currentTimeMillis());
					// delete permantently ReformatCBASE_Step1.main(new String[]{args[0],args[2],args[3]});
					//System.out.println(System.currentTimeMillis());
					ReformatCBASE.main(new String[]{args[0],args[1],args[2],args[3]});

					
					update_status("Estimate parameters for Bayesian model",8);
					System.out.println("Execute the parameter estimation for CBASE");
					//System.out.println(System.currentTimeMillis());
					CBASE_Solutions.main(new String[]{args[0],args[2]});
			
					String[] entities=entities(args[2]);
					//System.out.println(entities.length);
					//for (int i=0;i<entities.length;i++){
					//	System.out.println(entities[i]);
					//}
					//System.exit(0);
					
					boolean[] compute_uniform=new boolean[entities.length];
					ArrayList<String>[] samples=samples(args[2],entities);
					int[] count=count(args[0]+"AffinityCounts/AffinityCount.txt",samples);
					
					update_status("Computing significance",9);
					System.out.println("Computing Significance");
					//System.out.println(System.currentTimeMillis());
					for (int i=0;i<entities.length;i++){
						System.out.println(entities[i]);
						if(samples[i].size()<Integer.parseInt(args[6])||count[i]<Integer.parseInt(args[7])){
							compute_uniform[i]=true;
						}
						else{
							compute_uniform[i]=false;
						}
					}
	
					ComputeSignificance.main(new String[]{args[0],args[2],args[3]},entities,compute_uniform);
					
					
					
					
					
					//ArrayList<String>[] samples=samples(args[2],entities);
					//int[] count=count(args[0]+"AffinityCounts/AffinityCount.txt",samples);
					
					
					
					
					System.out.println("Start Filtering of significant genes");
					update_status("Significance Filtering Step 1",10);
					System.out.println("Filtering Step 1 - Preparing Blat Queries");
					//System.out.println(System.currentTimeMillis());
					Filter_Step1.main(new String[]{args[0],args[2],args[3]},entities);
	
					
					update_status("Significance Filtering Step 2",11);
					System.out.println("Performing Blat Queries");
					//System.out.println(System.currentTimeMillis());
					//File[] files_query=new File(args[0]+"PostSignFilter/Queries/").listFiles();
					
					//mkdir(args[0]+"PostSignFilter/OutputsBLAT");
					
					
					
					//String output_file_blat=args[0]+"PostSignFilter/OutputBLAT.txt";
					String command_blat="";//"cd "+args[3]+"SignificanceFilter/"+" && ./blat hg19.2bit "+args[0]+"PostSignFilter/Query.fa"+" -ooc=11.ooc -out=pslx "+args[0]+"PostSignFilter/OutputBLAT.txt";//+" && cd "+args[0];
				
					String OS = System.getProperty("os.name").toLowerCase();
					if(OS.indexOf("mac") >= 0){
						command_blat="cd "+args[3]+"SignificanceFilter/"+" && ./blat_mac hg19.2bit "+args[0]+"PostSignFilter/Query.fa"+" -ooc=11.ooc -out=pslx "+args[0]+"PostSignFilter/OutputBLAT.txt";//+" && cd "+args[0];	
					}
					else if(OS.indexOf("nix") >= 0 || OS.indexOf("nux") >= 0 || OS.indexOf("aix") > 0){
						command_blat="cd "+args[3]+"SignificanceFilter/"+" && ./blat_linux hg19.2bit "+args[0]+"PostSignFilter/Query.fa"+" -ooc=11.ooc -out=pslx "+args[0]+"PostSignFilter/OutputBLAT.txt";//+" && cd "+args[0];	
					}
					
					if(!command_blat.equals("")){
						execute(command_blat);
					}
					
					
					update_status("Significance Filtering Step 3",12);
					if(new File(args[0]+"/PostSignFilter/OutputBLAT.txt").exists()){
						System.out.println("Filtering Step 2 - Evaluating Blat Results");
						//System.out.println(System.currentTimeMillis());
						Filter_Step2.main(new String[]{args[0],args[2]},entities,compute_uniform);
						
					}
					
					
					update_status("Significance Filtering Step 4",13);
					System.out.println("Filtering Step 3 - Filtering out Genes");
					//System.out.println(System.currentTimeMillis());
					Filter_Step3.main(new String[]{args[0],args[2],args[3]},entities,compute_uniform);
					update_status("Completed",14);
					success=true;
				}
				catch(Exception e){
					System.out.println(e);
					success=false;
				}
				System.setOut(ps_console);
				//System.out.println(success);
				/*
				try{
					//System.out.println("START");
					
					
					
					
					update_status("Step1", 1);
					Thread.sleep(1000);
					update_status("Step2", 2);
					Thread.sleep(1000);
					update_status("Step3", 6);
					Thread.sleep(1000);
					update_status("Step4", 10);
					
					//System.out.println("DONE");
					success=true;
					ii=page_completion;
					reset();
					paint();
				}
				catch(Exception e){
					
				}*/
				is_running=false;
				ii=page_completion;
				reset();
				paint();
				
			}
		}
		
		public static void reindex_samples(String file_in, String file_out){
			try{
				FileInputStream in=new FileInputStream(file_in);
				DataInputStream inn=new DataInputStream(in);
				BufferedReader input= new BufferedReader(new InputStreamReader(inn));
				FileWriter out=new FileWriter(file_out);
				BufferedWriter output= new BufferedWriter(out);
				String[] tt=input.readLine().split("	");
				int[] index={index("Sample",tt),index("Cohort",tt)};
				output.write("ID	Sample	Cohort");
				output.newLine();
				int n=0;
				String s="";
				while((s=input.readLine())!=null){
					String[] t=s.split("	");
					output.write(n+"	"+t[index[0]]+"	"+t[index[1]]);
					output.newLine();
					n++;
				}
				input.close();
				output.close();
			}
			catch(Exception e){
				System.out.println(e);
			}
		}
		
		public static void execute(String command){
			System.out.println("starting "+command);//+"	"+system_cpu_load()+"	"+process_cpu_load());
			//boolean isWindows = System.getProperty("os.name").toLowerCase().startsWith("windows");
			
			//if(!new File(output_file).exists()){
				Process process=null;
				try{
					//if (isWindows) {
					//    process = Runtime.getRuntime().exec("cmd.exe /c "+command);
					//} else {
						String[] commands = { "/bin/bash", "-c", command };
						process = Runtime.getRuntime().exec(commands);
						
					//}	
					process.waitFor();
					
				}
				catch(Exception e){
					System.out.println(e);
				}
			//}
			
			System.out.println("done "+command);
		}
		
		public static ArrayList<String>[] samples(String file,String[] entities){
			ArrayList<String>[] samples=new ArrayList[entities.length];
			for (int i=0;i<samples.length;i++){
				samples[i]=new ArrayList<String>();
			}
			try{
				FileInputStream in=new FileInputStream(file);
				DataInputStream inn=new DataInputStream(in);
				BufferedReader input= new BufferedReader(new InputStreamReader(inn));
				int[] index_header=index_header(input.readLine().split("	"),index_header_samples);
				String s="";
				while((s=input.readLine())!=null){
					String entity=s.split("	")[index_header[2]];
					String sample=s.split("	")[index_header[1]];
					for (int i=0;i<entities.length;i++){
						if(entities[i].equals("PanCancer")||entities[i].equals(entity)){
							samples[i].add(sample);
						}
					}
				}
				input.close();
			}
			catch(Exception e){
				System.out.println(e);
			}
			return samples;
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
		
		public static String[] entities(String file){
			try{
				FileInputStream in=new FileInputStream(file);
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
				String [] entity=new String[aa.size()];
				for (int i=0;i<aa.size();i++){
					entity[i]=aa.get(i);
				}
				return entity;
			}
			catch(Exception e){
				System.out.println(e);
			}
			return new String[0];
		}
		
		public static int[] count(String file,ArrayList<String>[] samples){
			int[] counts=new int[samples.length];
			try{
				FileInputStream in=new FileInputStream(file);
				DataInputStream inn=new DataInputStream(in);
				BufferedReader input= new BufferedReader(new InputStreamReader(inn));
				input.readLine();
				String s="";
				while((s=input.readLine())!=null){
					String[] t=s.split("	");
					int count=Integer.parseInt(t[1])+Integer.parseInt(t[2])+Integer.parseInt(t[3])+Integer.parseInt(t[4])+Integer.parseInt(t[5])+Integer.parseInt(t[6]);
					for (int i=0;i<samples.length;i++){
						for (int j=0;j<samples[i].size();j++){
							if(samples[i].get(j).equals(t[0])){
								counts[i]+=count;
							}
						}
					}
				}
				input.close();
			}
			catch(Exception e){
				System.out.println(e);
			}
			return counts;
		}
		
		public static boolean contains (String s, ArrayList<String> t){
			for (int i=0;i<t.size();i++){
				if(t.get(i).equals(s)){
					return true;
				}
			}
			return false;
		}
		
		
		public void launch_mutpanning(){
			if(!path_hg19_file.substring(path_hg19_file.length()-1, path_hg19_file.length()).equals("/")&&!path_hg19_file.substring(path_hg19_file.length()-1, path_hg19_file.length()).equals("\\")){
				path_hg19_file=path_hg19_file+"/";
			}
			if(!path_temporary_file.substring(path_temporary_file.length()-1, path_temporary_file.length()).equals("/")&&!path_temporary_file.substring(path_temporary_file.length()-1, path_temporary_file.length()).equals("\\")){
				path_temporary_file=path_temporary_file+"/";
			}
			
			thread.start();
			
			
	        //return true;
		}
		
	}
	
	 public static void goWebsite(JLabel website, String url) {
	        website.addMouseListener(new MouseAdapter() {
	            public void mouseClicked(MouseEvent e) {
	            	//System.out.println(url);
	                try {
	                    Desktop.getDesktop().browse(new URI(url));
	                } catch (URISyntaxException | IOException ex) {
	                    //It looks like there's a problem
	                }
	            }
	        });
	   }
	
	 public static void goFolder(JLabel website, String folder) {
	        website.addMouseListener(new MouseAdapter() {
	            public void mouseClicked(MouseEvent e) {
	            	if (Desktop.isDesktopSupported()) {
	            	    try{
	            	    	Desktop.getDesktop().open(new File(folder));
	            	    }
	            	    catch(Exception ee){
	            	    	System.out.println(ee);
	            	    }
	            	}
	            }
	        });
	   }
	 
	private static class Bullet extends JComponent {
		boolean active=false;
		public Bullet(boolean active){
			this.active=active;
		}
		  public void paint(Graphics g) {
			  if(active){
				  g.setColor(Color.BLUE);
			  }
			  else{
				  g.setColor(Color.GRAY);
			  }
			  g.fillOval (0, 0, 10, 10);  
		  }
	}

	
	
}
