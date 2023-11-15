import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.apache.poi.hssf.usermodel.HSSFRow;
import org.apache.poi.hssf.usermodel.HSSFSheet;  
import org.apache.poi.hssf.usermodel.HSSFWorkbook;  
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.DataFormatter;
import org.apache.poi.ss.usermodel.FormulaEvaluator;  
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.util.CellReference;
import org.apache.poi.xssf.usermodel.XSSFRow;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;  
public class SortSimulation {
    public static List<Integer> BubbleSort(List<Integer> data) {
        for (int i = 0; i < data.size(); i++) {
            for (int j = 0; j < data.size(); j++) {
                int temp = data.get(j);
                if (j > 0 && temp < data.get(j - 1)) {
                    data.set(j, data.get(j - 1));
                    data.set(j - 1, temp);
                }
            }
        }
        return data;
    }

    public static List<Integer> InsertionSort(List<Integer> data) {
        for (int i = 1; i < data.size(); i++) {
        	int j = i;
        	int k = i - 1;
            while (k >= 0 && j >= 0) {
            	int temp = data.get(k);
            	if (temp > data.get(j)) {
		            data.set(k, data.get(j));
		            data.set(j, temp);
		            j = k;
	            }
            	k--;
            }
        }
        return data;
    }

    public static List<Integer> SelectionSort(List<Integer> data) {
        for (int i = 0; i < data.size(); i++) {
            int min = Integer.MAX_VALUE;
            int k = 0;
            for (int j = i; j < data.size(); j++) {
                if (min > data.get(j)) {
                    min = data.get(j);
                    k = j;
                }
            }
            data.set(k, data.get(i));
            data.set(i, min);
        }
        return data;
    }

    public static List<Integer> Merge(List<Integer> left, List<Integer> right) {
        List<Integer> data = new ArrayList<Integer>();
        int k = 0; // Index to represent right half.
        int i = 0; // Index to represent left half.
        while (i < left.size() || k < right.size()) {
            if (k < right.size() && i < left.size() && right.get(k) <= left.get(i)) {
                data.add(right.get(k++));
            } else if (i < left.size() && k < right.size() && left.get(i) < right.get(k)) {
                data.add(left.get(i++));
            } else if (i >= left.size()) {
                data.add(right.get(k++));
            } else if (k >= right.size()) {
                data.add(left.get(i++));
            }
        }
        return data;
    }

    public static List<Integer> MergeSort(List<Integer> data) {
        if (data.size() == 1) {
            return data;
        } 
        else {
            int m = Math.floorDiv(data.size() - 1, 2);
            List<Integer> L = data.subList(0,  m+1);
            List<Integer> R = data.subList(m+1, data.size());
            List<Integer> left = MergeSort(L);
            List<Integer> right = MergeSort(R);
            return Merge(left, right);
        }
    }
    
    /*public static List<List<Integer>> Partition(List<Integer> data) {
    	List<List<Integer>> lists = new ArrayList<List<Integer>>();
    	int pivot = Math.floorDiv(data.size()-1, 2);
    	int left_index = 0;
    	int right_index = data.size()-1;
    	int piv = data.get(pivot));
    	data.get(pivot, data.get(0));
    	data.get(0, piv);
    	while (right_index < mid_index) {
    		if (data.get(right_index) < data.get(0)) {
    			left_index++;
    		}
    		if (data.get(right_index) > data.get(0))
    			right_index--;
    		}
    		int temp = data.get(left_index);
    		data.set(left_index, data.get(right_index));
    		data.set(right_index, temp));
    		lists.add(data.subList(0, left_index+1));
    		lists.add(data.subList(right_index+1, data.size()));
    	}
    	return lists;
    }*/

   /* public static List<Integer> Partition(List<Integer> data, int start, int end) {
        int left_index = start;
        int piv_index = start;
        int right_index = end;
        List<Integer> temp = new ArrayList<Integer>();
        temp.addAll(data.subList(start, end));
        //System.out.println("SubList: " + temp);
        int pivot = Math.floorDiv(temp.size()-1, 2);
        for (int i = 0; i < temp.size(); i++) {
        	int value = temp.get(i);
        	if (value < temp.get(pivot)) {
        		if (piv_index > left_index /*&& piv_index + (piv_index-left_index) < right_index) {
        			for (int j = left_index; j < piv_index; j++) {
        				//Collections.swap(data, j, j+(piv_index-left_index));
        				if (j + piv_index-left_index < temp.size()) data.set(j+(piv_index-left_index), data.get(j));
        			}
        		}
        		data.set(left_index++, value);
        		piv_index++;
        	}
        	else if (value > temp.get(pivot)) {
        		data.set(--right_index, value);
        	}
        	else {
        		data.set(piv_index++, value);
        	}	
        }
        return new ArrayList<Integer>(Arrays.asList(left_index, piv_index, right_index));
    } */
    
    public static List<List<Integer>> Partition(List<Integer> data, int start, int end, int pivot) {
        int left_index = start;
        int piv_index = start;
        int right_index = end;
        List<Integer> temp = new ArrayList<Integer>();
        temp.addAll(data.subList(start, end));
        List<List<Integer>> lists = new ArrayList<List<Integer>>();
        for (int i = 0; i < temp.size(); i++) {
        	int value = temp.get(i);
        	if (value < pivot) {
        		if (piv_index > left_index) {
        			for (int j = left_index; j < piv_index; j++) {
        			if (j+(piv_index-left_index) < right_index) data.set(j+(piv_index-left_index), data.get(j));
        			}
        		}
        		data.set(left_index++, value);
            	piv_index++;
        	}
        	else if (value > pivot) {
        		data.set(--right_index, value);
        	}
        	else
        		data.set(piv_index++, value);
        }
        lists.add(data);
        lists.add(new ArrayList<Integer>(Arrays.asList(left_index)));
        lists.add(new ArrayList<Integer>(Arrays.asList(piv_index)));
        lists.add(new ArrayList<Integer>(Arrays.asList(right_index)));
        return lists;
    }
    
    public static List<Integer> MedianSort(List<Integer> data) {
    	List<Integer> medians = new ArrayList<Integer>();
    	int i = 0;
    	int n = data.size();
    	while (i < n) {
    		List<Integer> list;
    		if (i+11 < n) {
    			list = InsertionSort(data.subList(i, i+5));
    			int m = Math.floorDiv(list.size()-1, 2);
    			medians.add(list.get(m));
    			i += 11;
    		}
    		if (i+11 >= n) {
    			list = InsertionSort(data.subList(i, n));
    			int m = Math.floorDiv(list.size()-1, 2);
    			medians.add(list.get(m));
    			break;
    		}
    	}
    	return medians;
    }
    
    public static int MedianOfMedians(List<Integer> data) {
    	int m = Math.floorDiv(data.size()-1, 2);
    	return MedianOfMedians(data, 0, data.size(), m);
    }
    
    public static int MedianOfMedians(List<Integer> data, int start, int end, int k) {
    	List<Integer> sortedMedians = MedianSort(data);
    	int pivot = Math.floorDiv(sortedMedians.size()-1, 2);
    	List <List<Integer>> lists = Partition(sortedMedians, 0, sortedMedians.size(), sortedMedians.get(pivot));
    	sortedMedians = lists.get(0);
    	int left = lists.get(1).get(0);
    	int piv = lists.get(2).get(0);
    	int right = lists.get(3).get(0);
    	if (k == pivot) return sortedMedians.get(k);
    	else if (k < pivot) {
    		return MedianOfMedians(sortedMedians, 0, left, Math.floorDiv(left-1, 2));
    	}
    	else 
    		return MedianOfMedians(sortedMedians, right, sortedMedians.size(), Math.floorDiv(k-sortedMedians.size()-1+piv, 2));
    }
    
    public static int MedianHelper(List<Integer> data) {
    	return MedianOfMedians(data);
    }
     
    public static List<Integer> QuickSort(List<Integer> data, boolean median) {
    	return QuickSort(data, 0, data.size(), median);
    }
    
    public static List<Integer> QuickSort(List<Integer> data, int start, int end, boolean median) {
    	if (start == end) {
    		return data;
    	}
        else {
        	int pivot;
        	if (median) pivot = MedianHelper(data.subList(start, end));
        	else pivot = data.get(Math.floorDiv(start+end-1, 2));
        	List<List<Integer>> lists = Partition(data, start, end, pivot);
        	data = lists.get(0);
        	int left = lists.get(1).get(0);
        	int right = lists.get(3).get(0);
        	QuickSort(data, start, left, median);
            QuickSort(data, right, end, median);
            return data;	
        }
    } 
     
    public static double degreeOfSortedness(List<Integer> data) {
    	double inv = 0;
    	double n = data.size();
    	for (int j = 0; j < n; j++) {
    		for (int i = 0; i < j; i++) {
    			if (data.get(i) > data.get(j)) inv++;
    		}
    	}
    	double proba = (n * (n-1)) / 2;
    	return inv / proba;
    }
    
    public static List<Integer> adjustedSortedness(List<Integer> data, double threshold) {
    	Random randalf = new Random();
    	data = MergeSort(data);
    	Collections.reverse(data);
    	int m = (int)Math.ceil(data.size() * threshold);
    	int a = (int)Math.floor(randalf.nextDouble() * (1/2+threshold));
    	List<Integer> subList = MergeSort(data.subList(a, m));
    	for (int i = a; i < m; i++) {
    		data.set(i, subList.get(i-a));
    	}
    	return data;
    }
    
    public static List<Integer> getUniformSet(int seed, int size) {
    	List<Integer> uniform = new ArrayList<Integer>();
    	Random randalf = new Random();
    	for (int i = 0; i < size; i++) {
    		uniform.add((int)Math.floor(seed * randalf.nextDouble()));
    	}
    	return uniform;
    }
    
    public static List<Integer> getGaussianSet(int seed, int size) {
    	Random randalf = new Random();
    	List<Integer> gaussian = new ArrayList<Integer>();
    	for (int i = 0; i < size; i++) {
    		gaussian.add((int) Math.floor(seed * randalf.nextGaussian()));
    	}
    	return gaussian;		
    }

    public static void test(List<Integer> dataset) {
    	List<Integer> data = dataset;
        System.out.println("BubbleSort list: " + BubbleSort(data));
        data = dataset;
        System.out.println("InsertionSort list: " + InsertionSort(data));
        data = dataset;
        System.out.println("SelectionSort list: " + SelectionSort(data));
        data = dataset;
        System.out.println("MergeSort list: " + MergeSort(data));
        data = dataset;
        System.out.println("QuickSort list w/o median: " + QuickSort(data, false));
    }
    
    public static void calculateStats(List<Integer> dataset) {
    	System.gc();
    	List<Integer> data = dataset;
    	long startTime = System.nanoTime();
    	double startMemory = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory()) / 1024d / 1024d;
    	data = BubbleSort(data);
    	long endTime = System.nanoTime(); 
    	double runTime = (double) (endTime - startTime) / 1000000;
    	double endMemory = (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1024d / 1024d;
    	System.out.println("BubbleSort runtime: " + runTime + " milliseconds.");
    	System.out.println("BubbleSort memory usage: " + (endMemory - startMemory));
    	
    	data = dataset;
    	startTime = System.nanoTime();
    	startMemory = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory()) / 1024d / 1024d;
    	data = InsertionSort(data);
    	endTime = System.nanoTime();
    	runTime = (double)(endTime - startTime)/1000000;
    	endMemory = (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1024d / 1024d;
    	System.out.println("InsertionSort runtime: " + runTime + " milliseconds.");
    	System.out.println("InsertionSort memory usage: " + (endMemory - startMemory));
    	
    	data = dataset;
    	startTime = System.nanoTime();
    	startMemory = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory()) / 1024d / 1024d;
    	data = SelectionSort(data);
    	endTime = System.nanoTime();
    	runTime =  (double)(endTime - startTime)/1000000;
    	endMemory = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory()) / 1024d / 1024d;
    	System.out.println("SelectionSort runtime: " + runTime + " milliseconds.");
    	System.out.println("SelectionSort memory usage: " + (endMemory - startMemory));
    	
    	data = dataset;
    	startTime = System.nanoTime();
    	startMemory = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory()) / 1024d / 1024d;
    	data = MergeSort(data);
    	endTime = System.nanoTime();
    	endMemory = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory()) / 1024d / 1024d;
    	System.out.println("MergeSort runtime: " + runTime + " milliseconds.");
    	System.out.println("MergeSort memory usage: " + (endMemory - startMemory));
    	
    	data = dataset;
    	startTime = System.nanoTime();
    	startMemory = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory()) / 1024d / 1024d;
    	data = QuickSort(data, false);
    	endTime = System.nanoTime();
    	runTime = (double)(endTime - startTime)/1000000;
    	endMemory = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory()) / 1024d / 1024d;
    	System.out.println("QuickSort runtime w/o median: " + runTime + " milliseconds.");
    	System.out.println("QuickSort memory usage w/o median: " + (endMemory - startMemory));
    	
    	data = dataset;
    	startTime = System.nanoTime();
    	startMemory = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory()) / 1024d / 1024d;
    	data = QuickSort(data, true);
    	endTime = System.nanoTime();
    	runTime = (double)(endTime - startTime)/1000000;
    	endMemory = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory()) / 1024d / 1024d;
      	System.out.println(Runtime.getRuntime().totalMemory());
    	System.out.println(Runtime.getRuntime().freeMemory());
    	System.out.println("QuickSort runtime w/ median: " + runTime + " milliseconds.");
    	System.out.println("QuickSort memory usage w/ median: " + (endMemory - startMemory));
    }
    
    public static void uniformTestSize() {
    	System.out.println("Uniform size test: \n");
    	for (int i = 100; i <= 100000; i*=1.77827) {
    		System.out.println("Size of dataset: " + i);
    		calculateStats(getUniformSet(i, i));
    		System.out.println("\n");
    	}
    }
    
    public static void gaussianTestSize() {
    	System.out.println("Gaussian size test: \n");
    	for (int i = 100; i <= 100000; i*=(Math.pow(1.77827, 2))) {
    		System.out.println("Size of dataset: " + i);
    		calculateStats(getGaussianSet(i, i));
    		System.out.println("\n");
    	}
    }
    
    public static void dataTestSize(List<Integer> data, String file) throws IOException {
    	for (int i = 100; i <= data.size(); i*=Math.pow(1.42, 2)) {
    		System.out.println("Size of dataset: " + i);
    		data = getDataSet(file);
    		calculateStats(data.subList(0, i));
    		System.out.println("\n");
    	}
    }
    
    public static void dataTestSize2(List<Integer> data, String file) throws IOException {
    	for (int i = 1000; i <= data.size(); i*=2) {
    		System.out.println("Size of dataset: " + i);
    		data = getDataSet2(file);
    		calculateStats(data.subList(0, i));
    		System.out.println("\n");
    	}
    }
    
    public static void uniformTestDegree() {
    	System.out.println("Uniform sortedness test: \n");
    	double threshold = 0.05;
    	while (threshold <= 1) {
    		List<Integer> data = adjustedSortedness(getUniformSet(10000, 1000), threshold);
    		System.out.println("Degree of sortedness of dataset: " + degreeOfSortedness(data));
    		calculateStats(data);
    		threshold+=0.15;
    		System.out.println("Threshold: " + threshold);
    		System.out.println("\n");
    	}
    }
    
    public static void gaussianTestDegree() {
    	System.out.println("Gaussian sortedness test: \n");
    	double threshold = 0.05;
    	while (threshold <= 1) {
    		List<Integer> data = adjustedSortedness(getGaussianSet(10000, 1000), threshold);
    		System.out.println("Degree of sortedness of dataset: " + degreeOfSortedness(data));
    		calculateStats(data);
    		threshold+=0.15;
    		System.out.println("Threshold: " + threshold);
    		System.out.println("\n");
    	}
    }
    
    public static void dataTestDegree(String file) throws IOException {
    	double threshold = 0.05;
    	while (threshold <= 1) {
    		List<Integer> data = adjustedSortedness(getDataSet(file), threshold);
    		System.out.println("Degree of sortedness of dataset: " + degreeOfSortedness(data));
    		calculateStats(data);
    		threshold+=0.15;
    		System.out.println("Threshold: " + threshold);
    		System.out.println("\n");
    	}
    }
    
    public static void dataTestDegree2(String file) throws IOException {
    	double threshold = 0.05;
    	while (threshold <= 1) {
    		List<Integer> data = adjustedSortedness(getDataSet2(file).subList(0, 32000), threshold);
    		System.out.println("Degree of sortedness of dataset: " + degreeOfSortedness(data));
    		calculateStats(data.subList(0, 32000));
    		threshold+=0.15;
    		System.out.println("Threshold: " + threshold);
    		System.out.println("\n");
    	}
    }
    
    
    public static List<Integer> getDataSet(String file) throws IOException {
    	List<Integer> data = new ArrayList<Integer>();
    	FileInputStream fs = new FileInputStream(new File(file));
    	XSSFWorkbook wb = new XSSFWorkbook(fs);
    	XSSFSheet sheet = wb.getSheetAt(0);
    	//FormulaEvaluator formulaEvaluator=wb.getCreationHelper().createFormulaEvaluator();
    	for (int rowIndex = 5; rowIndex <= sheet.getLastRowNum(); rowIndex++) {
    		XSSFRow row = sheet.getRow(rowIndex);
    		Cell cell = null;
    		if (row != null) cell = row.getCell(CellReference.convertColStringToIndex("F"));
    		if (cell != null) data.add((int)cell.getNumericCellValue());
    	}
    	data.addAll(data);
    	return data;
    }
    
    public static List<Integer> getDataSet2(String file) throws IOException {
    	List<Integer> data = new ArrayList<Integer>();
    	FileInputStream fs = new FileInputStream(new File(file));
    	XSSFWorkbook wb = new XSSFWorkbook(fs);
    	XSSFSheet sheet = wb.getSheetAt(0);
    	DataFormatter formatter = new DataFormatter();
    	//FormulaEvaluator formulaEvaluator=wb.getCreationHelper().createFormulaEvaluator();
    	for (int rowIndex = 1; rowIndex <= sheet.getLastRowNum(); rowIndex++) {
    		XSSFRow row = sheet.getRow(rowIndex);
    		Cell cell = null;
    		if (row != null) cell = row.getCell(CellReference.convertColStringToIndex("S"));
    		if (cell != null) {
    			String val = formatter.formatCellValue(cell);
    			data.add((int)Integer.parseInt(val));
    		}
    	}
    	return data;
    }
    
    public static void officialTest() throws IOException {
    	//uniformTestSize();
    	//System.out.println("\n");
    	//gaussianTestSize();
    	//uniformTestDegree();
    	//gaussianTestDegree();
    	//List<Integer> data = getDataSet("C:\\economicdata1970-2018.xlsx"); Make sure file is in appropriate directory.
    	//List<Integer> data = getDataSet2("C:\\AP IB Assessments by Perf Lvl.xlsx");
    	//dataTestSize2(data, "C:\\AP IB Assessments by Perf Lvl.xlsx");
    	//data = getDataSet2("C:\\AP IB Assessments by Perf Lvl.xlsx");
    	//dataTestSize(data,"C:\\economicdata1970-2018.xlsx" );
    	//dataTestDegree("C:\\economicdata1970-2018.xlsx");
    	//dataTestDegree2("C:\\AP IB Assessments by Perf Lvl.xlsx");
    	//List<Integer> data = new ArrayList<Integer>(Arrays.asList(1, 12, 5, 26, 7, 14, 5, 7, 2, 10, -4, 7, 26, 50, 4, 11, 15, 13, 16));
    	//List<Integer> data = new ArrayList<Integer>(getGaussianSet(10000, 1000));
    	//System.out.println(QuickSort(data, true));
    	//System.out.println(MedianHelper(data));
    	//System.out.println(MergeSort(data).get(Math.floorDiv(data.size()-1, 2)));
    	//calculateStats(new ArrayList(Arrays.asList(1, 12, 5, 26, 7, 14, 5, 7, 2, 10, -4, 7, 26, 50, 4, 11, 15, 13, 16)));
    }
    
    public static void main(String[] args) throws IOException {
        officialTest();
    }
}