#include <iostream>
#include <math.h>
#include <stdlib.h>
#include<string.h>
#include<msclr\marshal_cppstd.h>
#include <ctime>// include this header 
#include<mpi.h>
#pragma once

#using <mscorlib.dll>
#using <System.dll>
#using <System.Drawing.dll>
#using <System.Windows.Forms.dll>
using namespace std;
using namespace msclr::interop;

int* inputImage(int* w, int* h, System::String^ imagePath) //put the size of image in w & h
{
	int* input;


	int OriginalImageWidth, OriginalImageHeight;

	//*********************************************************Read Image and save it to local arrayss*************************	
	//Read Image and save it to local arrayss

	System::Drawing::Bitmap BM(imagePath);

	OriginalImageWidth = BM.Width;
	OriginalImageHeight = BM.Height;
	*w = BM.Width;
	*h = BM.Height;
	int* Red = new int[BM.Height * BM.Width];
	int* Green = new int[BM.Height * BM.Width];
	int* Blue = new int[BM.Height * BM.Width];
	input = new int[BM.Height * BM.Width];
	for (int i = 0; i < BM.Height; i++)
	{
		for (int j = 0; j < BM.Width; j++)
		{
			System::Drawing::Color c = BM.GetPixel(j, i);

			Red[i * BM.Width + j] = c.R;
			Blue[i * BM.Width + j] = c.B;
			Green[i * BM.Width + j] = c.G;

			input[i * BM.Width + j] = ((c.R + c.B + c.G) / 3); //gray scale value equals the average of RGB values

		}

	}
	return input;
}


void createImage(int* image, int width, int height, int index)
{
	System::Drawing::Bitmap MyNewImage(width, height);


	for (int i = 0; i < MyNewImage.Height; i++)
	{
		for (int j = 0; j < MyNewImage.Width; j++)
		{
			//i * OriginalImageWidth + j
			if (image[i * width + j] < 0)
			{
				image[i * width + j] = 0;
			}
			if (image[i * width + j] > 255)
			{
				image[i * width + j] = 255;
			}
			System::Drawing::Color c = System::Drawing::Color::FromArgb(image[i * MyNewImage.Width + j], image[i * MyNewImage.Width + j], image[i * MyNewImage.Width + j]);
			MyNewImage.SetPixel(j, i, c);
		}
	}


	MyNewImage.Save("..//Data//Output//outputRes" + index + ".png");
	cout << "result Image Saved " << index << endl;
}


int main()
{
	int ImageWidth = 4, ImageHeight = 4;

	int start_s, stop_s, TotalTime = 0;

	System::String^ imagePath;
	std::string img;
	img = "..//Data//Input//test.png";

	imagePath = marshal_as<System::String^>(img);
	int* imageData = inputImage(&ImageWidth, &ImageHeight, imagePath); //pixel values of image

	start_s = clock();

	///start code here

	int kernel_size = 7;
	float** kernel = new float* [kernel_size];
	for (int i = 0; i < kernel_size; i++)
	{
		kernel[i] = new float[kernel_size];
		for (int j = 0; j < kernel_size; j++)
		{
			kernel[i][j] = 1.0 / (kernel_size * kernel_size * 1.0);
			//cout << kernel[i][j] << endl;
		}
	}
	//cout <<endl<< kernel[1][1] << endl;

	MPI_Init(NULL, NULL);///////////////////// mpiexec "HPC_ProjectTemplate.exe"
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int size_;
	MPI_Comm_size(MPI_COMM_WORLD, &size_);

	//start and end rows
	int sub_image_rows = ImageHeight / size_;
	int start = sub_image_rows * rank;
	int end_ = start + sub_image_rows;
	if (rank == size_ - 1)
		end_ = ImageHeight;
	//cout << "h " << ImageHeight << endl;
	//cout << "sub_image_rows  " << sub_image_rows << endl;
	cout << rank << ' ' << start << ' ' << end_ << ' ' << endl;

	int root_rank = 0;
	int local_rows = end_ - start;
	int local_size = (local_rows)*ImageWidth;
	int* local_imageData = new int[local_size];
	//kernel compution
	for (int i = start; i < end_; i++)//image rows
	{
		for (int j = 0; j < ImageWidth; j++)//image columns
		{
			//(i,j) * kernel
			int new_val = 0;
			int x = -(kernel_size / 2);
			int y = -(kernel_size / 2);
			for (int m = 0; m < kernel_size; m++)//kernel rows
			{
				for (int l = 0; l < kernel_size; l++)//kernel columns
				{

					if ((j + y >= ImageWidth) || (j + y < 0)) {//columns
						y++;
						continue;
					}
					if ((rank == 0 && (i + x < start)) || ((rank == size_ - 1 )&& (i -start+ x >= local_rows))) {
						y++;
						continue;
					}

					new_val += kernel[m][l] * imageData[(i + x) * ImageWidth + (j + y)];
					y++;
				}
				x++; y = -(kernel_size / 2);
			}
			
			local_imageData[(i - start) * ImageWidth + j] = new_val;///lesa
		}
	}


	if (rank == root_rank)
	{
		// Define the receive counts
		int* counts = new int[size_];
		int* displacement = new int[size_];
		for (int i = 0; i < size_; i++)
		{
			if (i == size_ - 1)
			{  
				counts[i] = (ImageWidth * (ImageHeight)) - (local_size * (size_-1));

				displacement[i] = (i) * local_size;
				
				break;
				
			}
			displacement[i] = i * local_size;
			counts[i] = local_size;
		}
		


		MPI_Gatherv(local_imageData, local_size, MPI_FLOAT, imageData, counts, displacement, MPI_FLOAT, root_rank, MPI_COMM_WORLD);
		
		
	}
	else
	{
		MPI_Gatherv(local_imageData, local_size, MPI_FLOAT, NULL, NULL, NULL, MPI_FLOAT, root_rank, MPI_COMM_WORLD);
	}
	if (rank == 0) {
		stop_s = clock();
		TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
		//createImage(imageData, ImageWidth, ImageHeight, 20);

		cout << "time: " << TotalTime << endl;

		createImage(imageData, ImageWidth, ImageHeight, 100);
	}
	
	
	MPI_Finalize();//////////////////////////////////////////////////////////////////
	///end
	

	free(imageData);
	return 0;

}



