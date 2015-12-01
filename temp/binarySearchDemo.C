#include<iostream>
 using namespace std;

void bubbleSort(int array[], int size){
      cout<<"  Input array is: "<<endl;
      for(int j=0; j<size; j++)
      {
       //Displaying Array
       cout<<"\t\t\tValue at "<<j<<" Index: "<<array[j]<<endl;
      }
      cout<<endl;
    // Bubble Sort Starts Here
     int temp;
     for(int i2=0; i2<size; i2++)
   {
     for(int j=0; j<size-1; j++)
     {
        //Swapping element in if statement
           if(array[j]>array[j+1])
       {
        temp=array[j];
        array[j]=array[j+1];
        array[j+1]=temp;
       }
     }
   }
   // Displaying Sorted array
      cout<<"  Sorted Array is: "<<endl;
     for(int i3=0; i3<size; i3++)
   {
    cout<<"\t\t\tValue at "<<i3<<" Index: "<<array[i3]<<endl;
   }
}// Sort Function Ends Here

int binarySearch(int array[],int size, int key){
         int start=1, end=size;
         int mid=(start+end)/2;

  while(start<=end&&array[mid]!=key){
        if(array[mid]<key){
          start=mid+1;
      }
     else{
          end=mid-1;
          }
       mid=(start+end)/2;
     }// While Loop End

   if(array[mid]==key)
    return mid; //Returnig to main
    else
   return -1;//Returnig to main

   cout<<"\n\n\n";
}// binarySearch Function Ends Here

void binarySearchDemo(){
      cout<<"Enter 5 numbers randomly : "<<endl;

      // Size can be change by replacing 5
      int array[5]; //Declaring array
     for(int i=0; i<5; i++)
      {
       cout<<"\t";  cin>>array[i]; // Initializing array
      }

      //Passing Arrary for Sorting
       bubbleSort(array,5);

    // Array has Sorted At This Point
    cout<<"\n\t\t\tEnter Key To Search: ";
    int key;
    cin>>key;

	//Passing Array, size and key To Search Key
	int result = binarySearch(array,5,key);
	cout << "Index of item: " << result << endl;
}
