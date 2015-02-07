#include "stdafx.h"
#include "Matrix.h"
#include <math.h>
#include <stdlib.h>
//#include "alldef.h"
#include <conio.h>
#include <stdio.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
CMatrix::CMatrix()
{
	m_nRow = 0;
	m_nCol = 0;

	m_pTMatrix.resize (m_nRow);
	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			m_pTMatrix[i].resize (m_nCol);
			m_pTMatrix[i][j] = (double) 0;
		}
	}
}

CMatrix::~CMatrix()
{	

}

CMatrix::CMatrix(unsigned int nRow,unsigned int nCol)
{
	// 动态分配二维数组
	TMatrix tMatrix;
	tMatrix.resize (nRow);

	for(unsigned int i=0; i < nRow; i++)
	{
		for(unsigned int j=0; j < nCol; j++)
		{
	        tMatrix[i].resize(nCol);
			tMatrix[i][j] = (long double) 0.00001;
		}
	}

	// 对对象变量赋值
	m_nRow	= nRow;
	m_nCol	= nCol;
	m_pTMatrix = tMatrix;

}

CMatrix::CMatrixone(unsigned int nRow,unsigned int nCol)
{
	// 动态分配二维数组
	TMatrix tMatrix;
	tMatrix.resize (nRow);

	for(unsigned int i=0; i < nRow; i++)
	{
		for(unsigned int j=0; j < nCol; j++)
		{
	        tMatrix[i].resize(nCol);
			tMatrix[i][j] = (long double) 0.001;
		}
	}

	// 对对象变量赋值
	m_nRow	= nRow;
	m_nCol	= nCol;
	m_pTMatrix = tMatrix;

}

CMatrix::CMatrix(CMatrix& cMatrixB)
{
	// Initialize the variable
	m_nRow = cMatrixB.m_nRow ;
	m_nCol = cMatrixB.m_nCol ;
	m_pTMatrix = cMatrixB.m_pTMatrix ;

	// Copy Data
	for(unsigned int i=0; i< cMatrixB.m_nRow; i++)
	{
		for(unsigned int j=0; j < cMatrixB.m_nCol; j++)
		{
			m_pTMatrix [i][j] = cMatrixB.m_pTMatrix [i][j];
		}
	}
	
}

CMatrix CMatrix::operator +(CMatrix& cMatrixB)
{
	// 要满足矩阵相加的条件: 行列数目相等!
	if(m_nRow != cMatrixB.m_nRow || m_nCol != cMatrixB.m_nCol )
	{
		::AfxMessageBox (TEXT("执行相加的两个矩阵维数不相等!"),MB_OK | MB_ICONERROR);
	}

	CMatrix	cMatrix = *this;

	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			cMatrix.m_pTMatrix [i][j] = m_pTMatrix [i][j] + cMatrixB.m_pTMatrix [i][j];
		}
	}

	return	cMatrix;

}

CMatrix CMatrix::operator -(CMatrix& cMatrixB)
{
	// 要满足矩阵相加的条件: 行列数数目相等!
	if(m_nRow != cMatrixB.m_nRow || m_nCol != cMatrixB.m_nCol )
	{
		::AfxMessageBox (TEXT("执行相减的两个矩阵维数不相等!"),MB_OK | MB_ICONERROR);
	}

	CMatrix cMatrix = *this;

	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			cMatrix.m_pTMatrix [i][j] = m_pTMatrix [i][j] - cMatrixB.m_pTMatrix [i][j];
		}
	}

	return	cMatrix;

}


CMatrix CMatrix::operator *(CMatrix& cMatrixB)
{
	if( m_nCol != cMatrixB.m_nRow )
	{
		::AfxMessageBox (TEXT("执行相乘的两个矩阵维数不满足相乘的条件!"),MB_OK | MB_ICONERROR);
	}
	
	CMatrix cResultMatrix(m_nRow,cMatrixB.m_nCol);

	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < cMatrixB.m_nCol; j++)
		{
			for(unsigned int m=0; m < m_nCol; m++)
			{
				cResultMatrix.m_pTMatrix [i][j] +=  m_pTMatrix [i][m] * cMatrixB.m_pTMatrix [m][j];
			}
		}
	}

	return cResultMatrix;
}


CMatrix CMatrix::operator * (double nValue)
{
	CMatrix cMatrix = *this;

	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			cMatrix.m_pTMatrix [i][j] =m_pTMatrix [i][j] * nValue;
		}
	}

	return cMatrix;
}


CMatrix& CMatrix::operator =(CMatrix& cMatrixB)
{
	if( (m_nRow != cMatrixB.m_nRow) || (m_nCol != cMatrixB.m_nCol) )
	{
		::AfxMessageBox(TEXT("等号左右两边的矩阵的维数不相等!"),MB_OK | MB_ICONERROR);
		return *this;	// return invalid value
	}	

	// 给变量赋值
	m_nRow = cMatrixB.m_nRow ;
	m_nCol = cMatrixB.m_nCol ;
	m_pTMatrix = cMatrixB.m_pTMatrix ;

	// 赋值操作
	for(unsigned int i=0; i < cMatrixB.m_nRow; i++)
	{
		for(unsigned int j=0; j< cMatrixB.m_nCol; j++)
		{
			m_pTMatrix [i][j] = cMatrixB.m_pTMatrix [i][j];
		}
	}
	
	return *this;
}

/*CMatrix& CMatrix::(CMatrix& cMatrixB)
{
	
	// 给变量赋值
	m_nRow = cMatrixB.m_nRow ;
	m_nCol = cMatrixB.m_nCol ;
	m_pTMatrix = cMatrixB.m_pTMatrix ;

	// 赋值操作
	for(unsigned int i=0; i < cMatrixB.m_nRow; i++)
	{
		for(unsigned int j=0; j< cMatrixB.m_nCol; j++)
		{
			m_pTMatrix [i][j] = cMatrixB.m_pTMatrix [i][j];
		}
	}
	
	return *this;
}*/

CMatrix& CMatrix::operator += (CMatrix& cMatrixB)
{
	if(m_nRow != cMatrixB.m_nRow || m_nCol != cMatrixB.m_nCol )
	{
		//printf("错误!执行相加的两个矩阵维数不相等!\n");
		::AfxMessageBox (TEXT("运算符的两边矩阵的维数不相等!"),MB_OK | MB_ICONERROR);
		return *this;	// return invalid value
	}
	
	// 赋值操作
	for(unsigned int i=0; i < cMatrixB.m_nRow; i++)
	{
		for(unsigned int j=0; j< cMatrixB.m_nCol; j++)
		{
			m_pTMatrix [i][j] += cMatrixB.m_pTMatrix [i][j];
		}
	}
	
	return *this;

}


CMatrix CMatrix::Transpose()
{
	CMatrix cMatrix(m_nCol,m_nRow);

	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			cMatrix.m_pTMatrix [j][i] = m_pTMatrix [i][j];
		}
	}

	return cMatrix;
}


CMatrix CMatrix::MergeColumnsToColumnVector()
{
	CMatrix cMatrix(m_nRow * m_nCol,(unsigned int)1);

	// 对矩阵赋值
	for(unsigned int j=0; j < m_nCol; j++)
	{
		for(unsigned int i=0; i < m_nRow; i++)
		{
			cMatrix.m_pTMatrix [i + j * m_nRow][(unsigned int)0] = m_pTMatrix [i][j];
		}
	}

	return cMatrix;

}


double CMatrix::GetTotalElementValue()
{
	double	nTotalValue = 0;

	for(unsigned int i=0; i < m_nRow; i++)
	{
		for( unsigned int j=0; j < m_nCol; j++)
		{
			nTotalValue += m_pTMatrix [i][j];
		}
	}

	return nTotalValue;
}

double	CMatrix::GetSystemError() const
{
	double	nSystemError = 0;

	for(unsigned int i=0; i < m_nRow; i++)
	{
		for( unsigned int j=0; j < m_nCol; j++)
		{
			nSystemError += m_pTMatrix [i][j] * m_pTMatrix [i][j];
		}
	}

	return nSystemError;

}


CMatrix CMatrix::AbsoluteValue ()
{
	CMatrix cMatrix = *this;

	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			cMatrix.m_pTMatrix [i][j] = fabs( m_pTMatrix [i][j]);

		}

	}
return cMatrix;
}


CMatrix CMatrix::Inverse()
{
	/////////////////////////////////////////////////////////////////////////
	// Using Gauss - Jordan Method
	// 参考书目: 计算机数值方法 --->施吉林 陈桂枝
	/////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////
	// 判断是否是可逆阵:
	//		可逆阵一定是方阵!!!

	if ( m_nRow != m_nCol)
	{
		//printf("错误!矩阵的行列数不相等,是非可逆阵!\n");
		::AfxMessageBox (TEXT("矩阵的行列数不相等,是非可逆阵!"),MB_OK | MB_ICONERROR);
	}

	CMatrix cMatrix = *this;

	//***********************************************************************
	// 思路:(非常规思维!)
	//		动态分配整型数组(2*m_nCol)来存储每次交换的行坐标的值
	//		不论有没有行交换都记录在数组中,
	//		1.没进行行交换的两个数据相等,在SwapMatrixRow()函数中
	//		检测到两个值相等立即返回,在SwapMatrixCol()函数中也一样,
	//		检测到两个值相等立即返回,不占用系统资源;
	//		2.不相等的就交换
	//***********************************************************************
    
	//	分配内存
	int *pIntArray = new int [2*m_nCol];

	// nSetp --- 约化步数,按列展开	
	for(unsigned int k=0; k < cMatrix.m_nCol; k++)
	{
		/////////////////////////////////////////////////////////////////////
		// 进行行交换 ---> 游戏规则:
		// 为保证计算过程的数值稳定性,在第k步约化时,先在{a(ik)}|i=k->n中选按
		// 模最大者作为约化主元素,并交换矩阵相应的行

		// 标记主元素
		double nMaxElement = cMatrix.m_pTMatrix [k][k];
		// 标记主元素所在的行数
		unsigned int nMainRow = k;

		for(unsigned int nCount = k+1; nCount < cMatrix.m_nCol; nCount++)
		{
			if( fabs(nMaxElement) < fabs(cMatrix.m_pTMatrix [nCount][k]) )
			{
				nMaxElement = cMatrix.m_pTMatrix [nCount][k];
				nMainRow = nCount;
			}
		}

		// 将欲交换的行数存在数组中
		pIntArray [2*k] = k;
		pIntArray [2*k+1] = nMainRow;
	

		// 交换行
		cMatrix.SwapMatrixRow(k,nMainRow);

		//Display();

		//	判断是否是可逆阵
		if(cMatrix.m_pTMatrix [k][k] == 0)
		{
			//printf("错误!此矩阵为非可逆阵!\n");
			::AfxMessageBox (TEXT("此矩阵为非可逆阵,没有逆矩阵!"),MB_OK | MB_ICONERROR);
		}

		cMatrix.m_pTMatrix [k][k] = 1/(cMatrix.m_pTMatrix [k][k]);
		

		// 算主列
		for(unsigned int i=0; i < cMatrix.m_nRow; i++)
		{	
			if( i != k)
				cMatrix.m_pTMatrix [i][k] = -(cMatrix.m_pTMatrix [k][k]) * (cMatrix.m_pTMatrix [i][k]); 
			
			//int nTempValue = m_pTMatrix [i][k];
			
		}
		
		//printf("\n");

		// 约化非主行
		for(unsigned int m=0; m < cMatrix.m_nRow; m++)
		{
			if ( m == k)
				continue;

			for(unsigned int n=0; n < cMatrix.m_nCol; n++)
			{
				if ( n == k)
					continue;

				cMatrix.m_pTMatrix [m][n] += cMatrix.m_pTMatrix [m][k] * cMatrix.m_pTMatrix [k][n];

				//printf("%10f ",m_pTMatrix [m][n]);

			}

			//printf("\n");

		}

		// 算主行
		for(unsigned int j=0; j < cMatrix.m_nCol; j++)
		{
			if( j != k)
				cMatrix.m_pTMatrix [k][j] = (cMatrix.m_pTMatrix [k][k]) * (cMatrix.m_pTMatrix [k][j]);

		}

	}

	
	/////////////////////////////////////////////////////////////////////////
	// 进行列交换 ---> 对交换行后的矩阵进行列交换 ---> 还原矩阵
	// 游戏规则:
	// 将开始矩阵中进行的行交换 ---> 现用相对应的列交换进行还原,即可得到所求的
	// 逆矩阵

	for(int i=2*m_nCol-1; i > 0; i--)
	{
		cMatrix.SwapMatrixCol(pIntArray[i],pIntArray[i-1]);
		i--;
	}

	delete []pIntArray;

	return cMatrix;

}

void CMatrix::SwapMatrixRow(unsigned int nRow1,unsigned int nRow2)
{
	if( nRow1 == nRow2)
		return;

	double *pArray = new double;

	for(unsigned int i=0; i < m_nCol; i++)
	{
		// Swap the datum of the two rows
		pArray[0] = m_pTMatrix [nRow1][i];
		m_pTMatrix [nRow1][i] = m_pTMatrix [nRow2][i];
		m_pTMatrix [nRow2][i] = pArray[0];
	}

	delete pArray;
}


void CMatrix::SwapMatrixCol(unsigned int nCol1,unsigned int nCol2)
{
	if( nCol1 == nCol2)
		return;

	double *pArray = new double;
	for(unsigned int i=0; i < m_nRow; i++)
	{
		// Swap the datum of the two columns
		pArray[0] = m_pTMatrix [i][nCol1];
		m_pTMatrix [i][nCol1] = m_pTMatrix [i][nCol2];
		m_pTMatrix [i][nCol2] = pArray[0];
	}
	
	delete pArray;
}

bool CMatrix::LoadDataFromFile(CString& strFileName)
{
	CStdioFile dataFile;
	LPCTSTR	lpszFileName = "";

	// CString convert to LPCTSTR
	strFileName.TrimLeft ();
	strFileName.TrimRight ();
	//strFileName.Format (lpszFileName);

	lpszFileName = (LPCTSTR)strFileName;

	if(!dataFile.Open (lpszFileName,CFile::modeRead | CFile::typeText))
	{
		::AfxMessageBox (TEXT("不能打开要读取数据的文件!"),MB_OK | MB_ICONERROR);
		dataFile.Close ();
		return FALSE;
	}

	// 用来存储提取文本文件中一行的数据
	CString	strData;				 
	// 用来记录文本文件中一共有多少行数据?
	unsigned int	nRow = 0;		

	/////////////////////////////////////////////////////////////////////////
	// Step 1: 得到文件的行列数目并根据文本文件的行列数目来设置对象(矩阵)的行
	// 列数目
	//

	while(dataFile.ReadString (strData) != FALSE)
	{
		++nRow;
	}

	// 根据文本文件的数据的行数设置对象(矩阵)的行数
	m_nRow = nRow;
	SetMatrixRowNumber(m_nRow);

	// 重新定位当前文件指针到文件开头
	dataFile.SeekToBegin ();

	dataFile.ReadString (strData);
	strData.TrimLeft ();
	strData.TrimRight ();

	TCHAR	SPACE_CHARACTER = ' ';
	// 用来记录文本文件中一行有多少列?
	unsigned int	nCol = 0;						
	
	// 空格符在字符串中的位置索引
	int nIndex = 0;

	do
	{
		nIndex = strData.Find (SPACE_CHARACTER);

		++nCol;

		// 提取字符串的子字符串,即提取一个double型实数数据
		//CString strDoubleNumber = strData.Left (nIndex);	

		// 将字符串转换为double型实数
		//double RealNumber = atof(strDoubleNumber);

		//int nTempNumber = strData.GetLength ();

		strData = strData.Right (strData.GetLength () - nIndex -1);

		// 去掉多余的空格
		strData.TrimLeft ();

		// Use for debugging
		//int nTempNum = strData.GetLength ();

	}while(nIndex != -1);

	// 根据文本文件的数据的列数设置对象(矩阵)的列数
	m_nCol = nCol;
	SetMatrixColNumber(m_nCol);

	// End of Getting the Rows and Cols of the Text File
	/////////////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////////////////////////////////
	// Step 2: 根据文本文件中的数据对矩阵赋值,并检测每行的列数是否和第一行的
	// 列数相等,不相等提示出错信息
	//

	// 重新定位当前文件指针到文件开头
	dataFile.SeekToBegin ();

	// 对矩阵中的元素装入文本文件的数据
	for(unsigned int i=0; i < m_nRow; i++)
	{
		dataFile.ReadString (strData);
		strData.TrimLeft ();
		strData.TrimRight ();

		// 验证每行的列数是否与第一行的列数相等
		unsigned int nVerifyColNum = 0;

		do
		{
			nIndex = strData.Find (SPACE_CHARACTER);

			++nVerifyColNum;

			if(nIndex != -1)
			{
				// 提取字符串的子字符串,即提取一个double型实数数据
				CString strDoubleNumber = strData.Left (nIndex);	

				// 将字符串转换为double型实数
				double RealNumber = atof(strDoubleNumber);

				m_pTMatrix [i][nVerifyColNum - 1] = RealNumber;
				//int nTempNumber = strData.GetLength ();

				strData = strData.Right (strData.GetLength () - nIndex -1);

				// 去掉多余的空格
				strData.TrimLeft ();

				// Using for debugging
				//double nReadNumber = m_pTMatrix [i][nVerifyColNum - 1];

				// Using for debugging
				//int nTempNum = strData.GetLength ();
			}
			else
			{
				double RealNumber = atof(strData);

				m_pTMatrix [i][nVerifyColNum - 1] = RealNumber;
			}

		}while(nIndex != -1);

		if(nVerifyColNum != m_nCol)
		{
			CString strRowNumber;
			strRowNumber.Format("%d",i + 1);

			CString strColNumber;
			strColNumber.Format("%d",m_nCol);

			CString strVerifyColNumber;
			strVerifyColNumber.Format("%d",nVerifyColNum);

			CString strMessage =  CString(TEXT("文本文件第")) + strRowNumber + CString(TEXT("行一共有")) + strVerifyColNumber + CString(TEXT("列,与第一行中的列数")) + strColNumber + CString(TEXT("不相等!"));

			LPCTSTR lpszText = "";
			lpszText = (LPCTSTR)strMessage;

			//strMessage.FormatMessage (lpszText);
			//::AfxMessageBox (lpszText,MB_OK);


			::AfxMessageBox (lpszText,MB_OK | MB_ICONERROR);
			dataFile.Close ();
			return false;
		}


	}


	dataFile.Close ();

	return TRUE;
}


bool CMatrix::SaveDataToFile (CString& strFileName)
{
	CStdioFile dataFile;
	LPCTSTR	lpszFileName = "";

	// CString convert to LPCTSTR
	strFileName.TrimLeft ();
	strFileName.TrimRight ();

	lpszFileName = (LPCTSTR)strFileName;

	if(!dataFile.Open (lpszFileName,CFile::modeCreate | CFile::modeNoTruncate | CFile::modeWrite | CFile::typeText))
	{
		::AfxMessageBox(TEXT("不能创建文件!"),MB_OK | MB_ICONERROR);
		dataFile.Close ();
		return false;
	}
	
	dataFile.SeekToEnd ();

	// 将对象(矩阵)中的数据写进文件
	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			 CString strRealNumber;

			 strRealNumber.Format ("%.16f  ", m_pTMatrix [i][j]);

			 // Using for debugging
			 //double nReadNumber = m_pTMatrix [i][j];

			 char *pBuffer = new char[strRealNumber.GetLength()];

			 memcpy(pBuffer,strRealNumber,strRealNumber.GetLength());

			 dataFile.Write (pBuffer,strRealNumber.GetLength ());

		}
		
		if( i != m_nRow - 1)
		{
			//char ReturnNewline[] = "\r\n";
			char ReturnNewline[] = "\n";
			
			dataFile.Write (ReturnNewline, (sizeof(ReturnNewline) - 1)/sizeof(char));
		}

	}


	dataFile.Close ();
	return true;
}

void CMatrix::Eye()
{
	// Verify whether the rows is equal to the columns or not
	if(m_nRow != m_nCol)
	{
		::AfxMessageBox (TEXT("此矩阵的行列数不相等!不能转变为单位阵!"),MB_OK | MB_ICONERROR);
		return;
	}
	
	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			if(i == j)
			{
				m_pTMatrix [i][j] =	1;
			}
			else
			{
				m_pTMatrix [i][j] =	0;
			}
		}

	}


}
void CMatrix::SetMatrixRowNumber(unsigned int nRow)
{
	m_nRow = nRow;

	m_pTMatrix.resize (m_nRow);
	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			m_pTMatrix[i].resize (m_nCol);
			m_pTMatrix[i][j] = (long double) 0;
		}
	}

}

void CMatrix::CopyMatrix(CMatrix& cMatrix)
{
	m_nRow	= cMatrix.m_nRow ;
	m_nCol	= cMatrix.m_nCol ;

	m_pTMatrix	= cMatrix.m_pTMatrix ;

	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			m_pTMatrix [i][j] = cMatrix.m_pTMatrix [i][j];
		}

	}

}


void CMatrix::GetMatrixData(CMatrix& cMatrix, unsigned int nIndex)
{
	if(m_nCol != 1)
	{
		::AfxMessageBox (TEXT("拷贝的矩阵不是列向量!"),MB_OK | MB_ICONERROR);
		return;
	}

	if((m_nRow - nIndex) < (cMatrix.m_nRow * cMatrix.m_nCol))
	{
		::AfxMessageBox (TEXT("拷贝矩阵的空间容量不足!"),MB_OK | MB_ICONERROR);
		return;
	}

	for(unsigned int i=0; i < cMatrix.m_nRow; i++)
	{
		for(unsigned int j=0; j < cMatrix.m_nCol; j++)
		{
			m_pTMatrix [nIndex + i * cMatrix.m_nCol + j][0] = cMatrix.m_pTMatrix [i][j];
		}

	}

}


void CMatrix::SetMatrixData(CMatrix& cMatrix, unsigned int nIndex)
{
	// Verify whether the colunm number is 1
	if(m_nCol != 1)
	{
		::AfxMessageBox (TEXT("本矩阵对象不是列向量,不满足条件!"),MB_OK | MB_ICONERROR);
		return;
	}

	// Verify whether the number of the object element is enough to be copyed

	if((m_nRow - nIndex) < (cMatrix.m_nRow * cMatrix.m_nCol))
	{
		::AfxMessageBox (TEXT("对象中的元素数量不足!"),MB_OK | MB_ICONERROR);
		return;
	}


	for(unsigned int i=0; i < cMatrix.m_nRow; i++)
	{
		for(unsigned int j=0; j < cMatrix.m_nCol; j++)
		{
			cMatrix.m_pTMatrix [i][j] = m_pTMatrix [nIndex + i * cMatrix.m_nCol + j][0];

			// Using for debugging
			//unsigned int nIndexNumber = nIndex + i * cMatrix.m_nRow + j;
			//double nData = cMatrix.m_pTMatrix [i][j];

		}
	}

}

void CMatrix::SetMatrixRowData(CMatrix& cMatrix, unsigned int nIndex, unsigned int nRow)
{
	// Verify whether the column number is 1
	if(m_nCol != 1)
	{
		::AfxMessageBox (TEXT("本矩阵对象不是列向量,不满足条件!"),MB_OK | MB_ICONERROR);
		return;
	}

	// Verify whether the number of the object element is enough to be copyed
	if((m_nRow - nIndex) < cMatrix.m_nCol )
	{
		::AfxMessageBox (TEXT("对象的元素数量不足!"),MB_OK | MB_ICONERROR);
		return;
	}

	for(unsigned int i=0; i < cMatrix.m_nCol; i++)
	{
		cMatrix.m_pTMatrix [nRow][i] = m_pTMatrix [nIndex + i][(unsigned int)0];
	}

}

void CMatrix::GetMatrixRowData(CMatrix& cMatrix, unsigned int nIndex, unsigned int nRow)
{
	if(m_nCol != 1)
	{
		::AfxMessageBox (TEXT("拷贝的矩阵不是列向量!"),MB_OK | MB_ICONERROR);
		return;
	}

	if((m_nRow - nIndex) < cMatrix.m_nCol)
	{
		::AfxMessageBox (TEXT("拷贝矩阵的空间容量不足!"),MB_OK | MB_ICONERROR);
		return;
	}
	
	for(unsigned int i=0; i < cMatrix.m_nCol; i++)
	{
		m_pTMatrix [nIndex + i][(unsigned int)0] = cMatrix.m_pTMatrix [nRow][i];
	}

}


void CMatrix::SetMatrixColNumber(unsigned int nCol)
{
	m_nCol = nCol;

	m_pTMatrix.resize (m_nRow);
	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			m_pTMatrix[i].resize (m_nCol);
			m_pTMatrix[i][j] = ( long double) 0;
		}
	}

}
/////////////////////////////////////////////////////////////////////////////
// 对矩阵中所有的元素进行一次非线性变换:
//		变换后的值y与变换前的值的关系是:
//			y = f(x) = 1 / (1 + exp(-x))	( 0 < f(x) < 1)
//
/////////////////////////////////////////////////////////////////////////////

CMatrix CMatrix::Sigmoid()
{
	CMatrix cMatrix = *this;

	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			cMatrix.m_pTMatrix [i][j] = 1 / (1 + exp(-m_pTMatrix [i][j]));
		}

	}
	
	return cMatrix;
}


/////////////////////////////////////////////////////////////////////////////
// 对矩阵中所有的元素进行一次非线性变换:
//		变换后的值y与变换前的值的关系是:
//			y = tanh(x) = (1 - exp(-x)) / (1 + exp(-x))
//					 =  1 - 2 * exp(-x) / (1 + exp(-x))	( -1 < f(x) < 1)
//
/////////////////////////////////////////////////////////////////////////////

CMatrix CMatrix::tanh ()
{
	CMatrix cMatrix = *this;

	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			cMatrix.m_pTMatrix [i][j] = 1 - (2 * exp(-m_pTMatrix [i][j])) / (1 + exp(-m_pTMatrix [i][j]));
		}

	}
	
	return cMatrix;
}

/////////////////////////////////////////////////////////////////////////////
// 对矩阵中所有的元素进行一次非线性变换:
//		变换后的值y与变换前的值的关系是:
//			y = Tansig(x) = 2 / (1 + exp(-2 * x)) -1
/////////////////////////////////////////////////////////////////////////////

CMatrix CMatrix::Tansig()
{
	CMatrix cMatrix = *this;

	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			cMatrix.m_pTMatrix [i][j] = 2 / (1 + exp(- 2 * m_pTMatrix [i][j])) - 1;
		}
	}
	
	return cMatrix;

}
/////////////////////////////////////////////////////////////////////////////
// 对矩阵中所有的元素进行一次非线性变换:
//		变换后的值y与变换前的值的关系是:
//			y = f'(x) = (1 / (1 + exp(-x)))'	( 0 < f(x) < 1)
//			  = exp(-x)/((1 + exp(-x))*(1 + exp(-x)))
/////////////////////////////////////////////////////////////////////////////

CMatrix CMatrix::SigmoidDerivative()
{
	CMatrix cMatrix = *this;

	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			cMatrix.m_pTMatrix [i][j] = exp(-m_pTMatrix [i][j]) / ((1 + exp(-m_pTMatrix [i][j])) * (1 + exp(-m_pTMatrix [i][j])));
		}

	}
	
	return cMatrix;
}


/////////////////////////////////////////////////////////////////////////////
// 对矩阵中所有的元素进行一次非线性变换:
//		变换后的值y与变换前的值的关系是:
//			y = tanh'(x) = ((1 - exp(-x)) / (1 + exp(-x)))'	( -1 < f(x) < 1)
//					 = 	2*exp(-x)/((1 + exp(-x))*(1 + exp(-x)))
/////////////////////////////////////////////////////////////////////////////

CMatrix CMatrix::tanhDerivative()
{
	CMatrix cMatrix = *this;

	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			cMatrix.m_pTMatrix [i][j] = 2 * exp(-m_pTMatrix [i][j]) / ((1 + exp(-m_pTMatrix [i][j])) * (1 + exp(-m_pTMatrix [i][j])));
		}

	}

	return cMatrix;
}
/////////////////////////////////////////////////////////////////////////////
// 对矩阵中所有的元素进行一次非线性变换:
//		变换后的值y与变换前的值的关系是:
//			y = Tansig'(x) = (2 / (1 + exp(-2 * x)) -1)'
//				= (2 / (1 + exp(-2 * x)) -1) * (2 / (1 + exp(-2 * x)) -1) -1
/////////////////////////////////////////////////////////////////////////////

CMatrix CMatrix::TansigDerivative()
{
	CMatrix cMatrix = *this;

	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			cMatrix.m_pTMatrix [i][j] = (2 / (1 + exp(- 2 * m_pTMatrix [i][j])) - 1) * (2 / (1 + exp(- 2 * m_pTMatrix [i][j])) - 1) - 1;
		}
	}

	return cMatrix;

}

void CMatrix::SetMatrixRowAndCol(unsigned int nRow,unsigned int nCol)
{
	m_nRow = nRow;
	m_nCol = nCol;

	// 分配内存
	m_pTMatrix.resize (m_nRow);
	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			m_pTMatrix[i].resize (m_nCol);
			m_pTMatrix[i][j] = (long double) 0;
		}
	}

}

void CMatrix::Initialize()
{
	m_nRow = 0;
	m_nCol = 0;

	m_pTMatrix.resize (m_nRow);
	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			m_pTMatrix[i].resize (m_nCol);
			m_pTMatrix[i][j] = (long double) 0;
		}
	}

}

void CMatrix::InitializeZero()
{
	m_nRow = 0;
	m_nCol = 0;

	m_pTMatrix.resize (m_nRow);
	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			m_pTMatrix[i].resize (m_nCol);
			m_pTMatrix[i][j] = (long double) 0;
		}
	}

}
void CMatrix::nncpy(CMatrix& cMatrix,unsigned int nTimes)
{
	m_nRow = cMatrix.m_nRow ;
	m_nCol = cMatrix.m_nCol * nTimes;

	// 根据空间分配内存
	m_pTMatrix.resize (m_nRow);
	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			m_pTMatrix[i].resize (m_nCol);
			m_pTMatrix[i][j] = (double) 0;
		}
	}

	// 对矩阵赋值
	for(i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < nTimes; j++)
		{
			for(unsigned int k=0; k < cMatrix.m_nCol; k++)
			{
				m_pTMatrix [i][j * cMatrix.m_nCol + k] = cMatrix.m_pTMatrix [i][k];
			}
		}
	}

}
void CMatrix::RandomInitialize ()
{
	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			m_pTMatrix [i][j] = (long double)(rand() - (0.5*RAND_MAX)) / (0.5*RAND_MAX);
		}
	}
	
}

void CMatrix::CopySubMatrix(CMatrix& cMatrix,unsigned int nStartX,unsigned int nStartY)
{
	if((m_nRow  < cMatrix.m_nRow + nStartY ) | (m_nCol  < cMatrix.m_nCol + nStartX))
	{
		::AfxMessageBox (TEXT("被拷贝的矩阵维数小于要拷贝的矩阵所需要的维数!"),MB_OK | MB_ICONERROR);
		return;
	}

	for(unsigned int i=0;  i < cMatrix.m_nRow; i++)
	{
		for(unsigned int j=0; j < cMatrix.m_nCol; j++)
		{
			cMatrix.m_pTMatrix [i][j] = m_pTMatrix [nStartY + i][nStartX + j];	
		}
	}

}

void CMatrix::CopySubMatrixFromVector(CMatrix& cMatrix,unsigned int nIndex)
{
	if(m_nCol != 1)
	{
		::AfxMessageBox (TEXT("被拷贝的矩阵不是列向量!!!"),MB_OK | MB_ICONERROR);
		return;
	}

	for(unsigned int j=0; j < cMatrix.m_nCol; j++)
	{
		for(unsigned int i=0; i < cMatrix.m_nRow; i++)
		{
			cMatrix.m_pTMatrix [i][j] = m_pTMatrix [nIndex + j * cMatrix.m_nRow + i ][(unsigned int)0];
		}
	}

}

void CMatrix::nncpyi(CMatrix &cMatrix, unsigned int nTimes)
{
	m_nRow	=	cMatrix.m_nRow ;
	m_nCol	=	cMatrix.m_nCol *	nTimes;

	// 根据空间分配内存
	m_pTMatrix.resize (m_nRow);
	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			m_pTMatrix[i].resize (m_nCol);
			m_pTMatrix[i][j] = (BYTE) 0;
		}
	}

	// 赋值
	for(i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < cMatrix.m_nCol; j++)
		{
			for(unsigned int k=0; k < nTimes; k++)
			{
				m_pTMatrix [i][j * nTimes + k] = cMatrix.m_pTMatrix [i][j];
			}
		}
	}

}

void CMatrix::nncpyd(CMatrix &cMatrix)
{
	m_nRow	=	cMatrix.m_nRow ;
	m_nCol	=	cMatrix.m_nCol * cMatrix.m_nRow ;

	// 根据空间分配内存
	m_pTMatrix.resize (m_nRow);
	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			m_pTMatrix[i].resize (m_nCol);
			m_pTMatrix[i][j] = (BYTE) 0;
		}
	}

	// 给矩阵赋值
	for(i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < cMatrix.m_nCol; j++)
		{
			for(unsigned int k=0; k < cMatrix.m_nRow; k++)
			{
				if(i == (j * cMatrix.m_nRow + k) % cMatrix.m_nRow )
					m_pTMatrix [i][j * cMatrix.m_nRow + k] = cMatrix.m_pTMatrix [i][j];
			}
		}
	}

}




void CMatrix::MakeAllColumnElementsSameValue(unsigned int nRowIndex)
{
	for(unsigned int i=0; i < m_nRow; i++)
	{
		if(i == nRowIndex)
			continue;

		for(unsigned int j=0; j < m_nCol; j++)
		{
			m_pTMatrix [i][j] = m_pTMatrix [nRowIndex][j];
		}
	}

}

CMatrix CMatrix::operator / (CMatrix& cMatrixB)
{
	CMatrix cMatrix = *this;

	if( (m_nRow != cMatrixB.m_nRow) || (m_nCol != cMatrixB.m_nCol) )
	{
		::AfxMessageBox (TEXT("两个矩阵的维数不相等,不满足矩阵点乘的条件!"),MB_OK | MB_ICONERROR);
		return cMatrix;	// return a invalid value
	}
	
	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			cMatrix.m_pTMatrix [i][j] = m_pTMatrix [i][j] * cMatrixB.m_pTMatrix [i][j];
		}

	}

	return cMatrix;

}

CMatrix WINAPI operator - (double nValue,CMatrix& cMatrixB)
{
	CMatrix	cMatrix = cMatrixB;

	for(unsigned int i=0; i < cMatrix.GetMatrixRowNumber (); i++)
	{
		for(unsigned int j=0; j < cMatrix.GetMatrixColNumber (); j++)
		{
			cMatrix.m_pTMatrix [i][j] = nValue - cMatrixB.m_pTMatrix [i][j];
		}
	}

	return cMatrix;
}

CMatrix WINAPI MergeMatrix(CMatrix& cMatrixA,CMatrix& cMatrixB)
{
	//	条件检测
	if( cMatrixA.GetMatrixRowNumber () != cMatrixB.GetMatrixRowNumber () )
	{
		::AfxMessageBox (TEXT("参与合并的两个矩阵的行数不相等!"),MB_OK | MB_ICONERROR);
	
		return cMatrixA;	// return invalid value
	}
	
		//CMatrix cMatrix = *this;

	CMatrix cMatrix(cMatrixA.GetMatrixRowNumber (),cMatrixA.GetMatrixColNumber () + cMatrixB.GetMatrixColNumber ());

	for(unsigned int i=0; i < cMatrixA.GetMatrixRowNumber (); i++)
	{
		for(unsigned int j=0; j < cMatrixA.GetMatrixColNumber (); j++)
		{
			cMatrix.m_pTMatrix [i][j] = cMatrixA.m_pTMatrix [i][j];
		}

		for(unsigned int k=0; k < cMatrixB.GetMatrixColNumber (); k++)
		{
			cMatrix.m_pTMatrix [i][cMatrixA.GetMatrixColNumber () + k] = cMatrixB.m_pTMatrix [i][k];
		}

	}


	return cMatrix;
}

void WINAPI LMForwardCalculateInit( int nInputLayerNumber,
													int nHideLayerNumber,
													int nOutputLayerNumber,
													CMatrix &matrixDemoDataInput,
													CMatrix &matrixInputLayerValue,
													CMatrix &matrixInputToHideWeightValue,
													CMatrix &matrixHideLayerValveValue,
													CMatrix &matrixHideToOutputWeightValue,
													CMatrix &matrixOutputLayerValveValue
												   )
{

	/************************************************************************
	*					--------->Use Matlab Method <---------              *
	************************************************************************/

		
	/////////////////////////////////////////////////////////////////////////
	// 构造输入层元素的矩阵
	// 构造规则:
	//		1. 样本数目做为矩阵的行数;
	//		2. 单个样本的输入层的数目做为矩阵的列数;
	//		3. 矩阵中的元素即为对应的输入层的值
	//

	CMatrix cMatrixInputLayerValue(matrixDemoDataInput.GetMatrixRowNumber (), nInputLayerNumber);

	// 得到样本的输入值
	matrixDemoDataInput.CopySubMatrix (cMatrixInputLayerValue,(unsigned int)0,(unsigned int)0);

	CMatrix cMatrixTInputLayerValue = cMatrixInputLayerValue.Transpose ();

	matrixInputLayerValue.CopyMatrix (cMatrixTInputLayerValue);


	/////////////////////////////////////////////////////////////////////////
	// 构造权值矩阵 -----> 由单个样本输入层与隐含层之间的权值做为元素组成
	// 构造规则:
	//		1. 单个样本的输入层的数目做为矩阵行数;
	//		2. 单个样本的隐含层的数目做为矩阵的列数;
	//		3. 对矩阵中的元素进行随机初始化,值在(-1,1)之间;
	//		4. 所有样本的输入层和隐含层的数目是相等的;
	//		5. 所有样本使用的是同一个权值矩阵.
	//

	CMatrix	cMatrixInputToHideWeightValue(nHideLayerNumber, nInputLayerNumber);

	// 随机初始化矩阵内元素的值
	cMatrixInputToHideWeightValue.RandomInitialize ();

	matrixInputToHideWeightValue.CopyMatrix (cMatrixInputToHideWeightValue);


	/////////////////////////////////////////////////////////////////////
	//	构造样本隐含层的阀值矩阵
	//	构造规则:
	//		1. 样本的数目做为矩阵行数;
	//		2. 单个样本的隐含层的数目做为矩阵的列数;
	//		3. 随机初始化矩阵中的元素,使值在(-1,1)之间;
	//		4. 矩阵中每行的数据都和第一行数据相对应的位置相等.
	//

	CMatrix cMatrixHideLayerValveValue(nHideLayerNumber,(unsigned int)1);

	// 随机初始化矩阵内元素的值
	cMatrixHideLayerValveValue.RandomInitialize ();

	matrixHideLayerValveValue.CopyMatrix(cMatrixHideLayerValveValue);


	/////////////////////////////////////////////////////////////////////
	//	构造权值矩阵	----->	由单个样本的隐含层与输出层之间权值做为元素
	//							组成
	//	构造规则:
	//		1. 单个样本的隐含层的数目做为矩阵的行数;
	//		2. 单个样本的输出层的数目做为矩阵的列数;
	//		3. 对矩阵中的元素进行随机初始化,值在(-1,1)之间;
	//		4. 所有样本的隐含层和输出层的数目是相等的;
	//		5. 所有样本使用的是同一个权值矩阵.
	//

	CMatrix cMatrixHideToOutputWeightValue(nOutputLayerNumber, nHideLayerNumber);

	// 对矩阵的元素随机初始化
	cMatrixHideToOutputWeightValue.RandomInitialize ();

	matrixHideToOutputWeightValue.CopyMatrix (cMatrixHideToOutputWeightValue);


	/////////////////////////////////////////////////////////////////////
	//	构造样本的输出层的阀值矩阵
	//	构造规则:
	//		1. 样本的数目做为矩阵的行数;
	//		2. 单个样本的输出层的数目做为矩阵的列数;
	//		3. 随机初始化矩阵中的元素,使值在(-1,1)之间;
	//		4. 矩阵中每行的数据都和第一行数据相对应的位置相等.
	//

	CMatrix cMatrixOutputLayerValveValue(nOutputLayerNumber,(unsigned int)1);

	// 随机初始化矩阵内元素的值
	cMatrixOutputLayerValveValue.RandomInitialize ();

	matrixOutputLayerValveValue.CopyMatrix(cMatrixOutputLayerValveValue);


}
/////////////////////////////////////////////////////////////////////////////
//	Levenberg-Marquart		---->		前向计算						   //
/////////////////////////////////////////////////////////////////////////////
 void WINAPI LMForwardCalculate ( int nInputLayerNumber,
												int nHideLayerNumber,
												int nOutputLayerNumber,
												//bool bSimulateDataFlag,
												//int nComboFunc,
												CMatrix &matrixDemoDataInput,
												CMatrix &matrixInputLayerValue,
												CMatrix &matrixInputToHideWeightValue,
												CMatrix &matrixHideLayerValveValue,
												CMatrix &matrixHideLayerOutput,
												CMatrix &matrixHideToOutputWeightValue,
												CMatrix &matrixOutputLayerOutput,
												CMatrix &matrixOutputLayerValveValue
											   )
{

	
		CMatrix cMatrixInputLayerValue(matrixDemoDataInput.GetMatrixRowNumber (), nInputLayerNumber);

		// 得到样本的输入值
		matrixDemoDataInput.CopySubMatrix (cMatrixInputLayerValue, (unsigned int)0, (unsigned int)0);

		CMatrix cMatrixTInputLayerValue = cMatrixInputLayerValue.Transpose ();

//cMatrixTInputLayerValue =cMatrixTInputLayerValue.Tansig();
		matrixInputLayerValue.CopyMatrix (cMatrixTInputLayerValue);





	/////////////////////////////////////////////////////////////////////////
	// 得到所有样本的隐含层的净输入
	// 构造规则:
	//		1. 样本的数目做为矩阵行数;
	//		2. 单个样本的隐含层的数目做为矩阵的列数;
	//		3. 矩阵元素中的值即为对应的样本的隐含层的净输入:
	//		   由 
	//				cMatrixInputLayerValue * cMatrixInputToHideWeightValue 
	//				+ cMatrixHideLayerValveValue
	//			得到.
	//
		
	CMatrix cMatrixExHideLayerValveValue;
	cMatrixExHideLayerValveValue.nncpyi (matrixHideLayerValveValue, matrixDemoDataInput.GetMatrixRowNumber ());

	CMatrix cMatrixHideLayerPureInput(nHideLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());

	cMatrixHideLayerPureInput = matrixInputToHideWeightValue * matrixInputLayerValue;

	cMatrixHideLayerPureInput += cMatrixExHideLayerValveValue;


	/////////////////////////////////////////////////////////////////////
	//	算出所有样本的隐含层的输出
	//	构造规则:
	//		1. 隐含层的输出y与隐含层的输入x的关系可用函数表示
	//			y = f(x)
	//		2. 矩阵的维数和隐含层的净输入矩阵的维数相等
	//

	CMatrix cMatrixHideLayerOutput(nHideLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());

	
		cMatrixHideLayerOutput = cMatrixHideLayerPureInput.Tansig();
	

	
	matrixHideLayerOutput.CopyMatrix(cMatrixHideLayerOutput);

	/////////////////////////////////////////////////////////////////////
	//	得到所有样本输出层的净输入
	//	构造规则;
	//		1. 样本的数目做为矩阵的行数;
	//		2. 单个样本的输出层的数目做为矩阵的列数;
	//		3. 矩阵中元素的值即为对应样本的输出层的净输入:
	//			由
	//				cMatrixHideLayerOutput * cMatrixHideToOutputWeightValue
	//				+ cMatrixOutputLayerValveValue
	//			得到
	//

	CMatrix cMatrixExOutputLayerValveValue;
	cMatrixExOutputLayerValveValue.nncpyi (matrixOutputLayerValveValue, matrixDemoDataInput.GetMatrixRowNumber ());

	CMatrix cMatrixOutputLayerPureInput(nOutputLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());

	cMatrixOutputLayerPureInput = matrixHideToOutputWeightValue * cMatrixHideLayerOutput;
	cMatrixOutputLayerPureInput  += cMatrixExOutputLayerValveValue;


	/////////////////////////////////////////////////////////////////////
	//	算出所有样本的输出层的输出
	//	构造规则:
	//		1. 矩阵的维数与得到的所有样本的输出层的净输入组成的矩阵一样;
	//		2. 输出层的输出y和输出层的输入可用关系式
	//			y = f(x)
	//		表示
	//

	CMatrix cMatrixOutputLayerOutput(nOutputLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());

	cMatrixOutputLayerOutput = cMatrixOutputLayerPureInput.Tansig ();
	



	matrixOutputLayerOutput.CopyMatrix(cMatrixOutputLayerOutput);

}
 bool WINAPI LMDemoDataTrainRepeat (	int nInputLayerNumber,
								int nHideLayerNumber,
								int nOutputLayerNumber,
								double nSystemErrorOld,
								double nSystemErrorNew,
								double nSystemErrorLevel,
								double nSystemError,
								double nStep,
								UINT nMaxTrainTimes,
								UINT nTrainTimes,
								CMatrix &matrixDemoDataInput,
								CMatrix &matrixInputLayerValue,
								CMatrix &matrixInputToHideWeightValue,
								CMatrix &matrixHideLayerValveValue,
								CMatrix &matrixHideLayerOutput,
								CMatrix &matrixHideToOutputWeightValue,
								CMatrix &matrixOutputLayerOutput,
								CMatrix &matrixOutputLayerValveValue
								)
{

	 ::LMForwardCalculate (nInputLayerNumber,
						nHideLayerNumber,
						nOutputLayerNumber,
						matrixDemoDataInput,
						matrixInputLayerValue,
						matrixInputToHideWeightValue,
						matrixHideLayerValveValue,
						matrixHideLayerOutput,
						matrixHideToOutputWeightValue,
						matrixOutputLayerOutput,
						matrixOutputLayerValveValue
						);
		

	/////////////////////////////////////////////////////////////////////
	//	算出所有样本的输出层的delta矩阵
	//	构造规则:
	//		1. 样本的数目为矩阵的行数;
	//		2. 样本输出层的数目为矩阵的列数;
	//		3. 矩阵中的元素的值y为:
	//				y = -(前向计算出的输出层的值 - 样本的输出层的值) * f'(net)
	//


	CMatrix cMatrixTDemoOutput(matrixDemoDataInput.GetMatrixRowNumber (), nOutputLayerNumber);
		
	// 得到样本中输出层的数据????
	matrixDemoDataInput.CopySubMatrix (cMatrixTDemoOutput, (unsigned int)nInputLayerNumber, (unsigned int)0);

	CMatrix cMatrixDemoOutput = cMatrixTDemoOutput.Transpose ();

	// 得到样本中输出层的误差
	CMatrix cMatrixOutputLayerError(nOutputLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());

	cMatrixOutputLayerError = cMatrixDemoOutput - matrixOutputLayerOutput;

		
	nSystemErrorOld = cMatrixOutputLayerError.GetSystemError ();

	for(int nLoopTimes=1; nLoopTimes < nMaxTrainTimes; nLoopTimes++)	
	{
		if(nSystemErrorOld < nSystemErrorLevel)
		{
			nLoopTimes--;
			break;
		}

		CMatrix	cMatrixExHideLayerOutput;
		cMatrixExHideLayerOutput.nncpyi (matrixHideLayerOutput, nOutputLayerNumber);
//6x22-->6x66

		CMatrix	cMatrixOutputLayerDelta (nOutputLayerNumber, matrixDemoDataInput.GetMatrixRowNumber());
		//3x22
		// 注意: 此处'/' 是 '点乘'!!!
		cMatrixOutputLayerDelta = 1 - matrixOutputLayerOutput / matrixOutputLayerOutput;
//3x22
		CMatrix cMatrixExOutputLayerDelta;
		cMatrixExOutputLayerDelta.nncpyd (cMatrixOutputLayerDelta);
		//3x66
		cMatrixExOutputLayerDelta = cMatrixExOutputLayerDelta * (-1.0);

	
		CMatrix cMatrixTHideToOutputWeightValue (matrixHideToOutputWeightValue.GetMatrixColNumber(), matrixHideToOutputWeightValue.GetMatrixRowNumber());
		cMatrixTHideToOutputWeightValue = matrixHideToOutputWeightValue.Transpose();
//6x3
		CMatrix cMatrixExHideLayerDelta;
		// 注意: 此处'/' 是 '点乘'!!!
		cMatrixExHideLayerDelta.CopyMatrix ( (1 - (cMatrixExHideLayerOutput / cMatrixExHideLayerOutput)) / ( cMatrixTHideToOutputWeightValue * cMatrixExOutputLayerDelta) );
			//???	

		CMatrix	cMatrixExInputLayerValue;
		cMatrixExInputLayerValue.nncpyi (matrixInputLayerValue, nOutputLayerNumber);

		
		CMatrix	cMatrixJ11;
		cMatrixJ11.nncpy (cMatrixExHideLayerDelta.Transpose(), cMatrixExInputLayerValue.GetMatrixRowNumber ());

		CMatrix cMatrixJ12;
		cMatrixJ12.nncpyi(cMatrixExInputLayerValue.Transpose (), cMatrixExHideLayerDelta.GetMatrixRowNumber());


		CMatrix	cMatrixJ1;
		// 注意: 此处'/' 是 '点乘'!!!
		cMatrixJ1.CopyMatrix (cMatrixJ11 / cMatrixJ12);


		CMatrix cMatrixJ21;
		cMatrixJ21.nncpy (cMatrixExOutputLayerDelta.Transpose (), cMatrixExHideLayerOutput.GetMatrixRowNumber ());

		CMatrix cMatrixJ22;
		cMatrixJ22.nncpyi (cMatrixExHideLayerOutput.Transpose (), cMatrixExOutputLayerDelta.GetMatrixRowNumber ());

		CMatrix cMatrixJ2;
		// 注意: 此处'/' 是 '点乘'!!!
		cMatrixJ2.CopyMatrix (cMatrixJ21 / cMatrixJ22);

		
		CMatrix cMatrixZ;
		cMatrixZ.CopyMatrix (MergeMatrix(MergeMatrix(MergeMatrix(cMatrixJ1, cMatrixExHideLayerDelta.Transpose ()), cMatrixJ2), cMatrixExOutputLayerDelta.Transpose ()));


		CMatrix cMatrixMOutputLayerError;
		cMatrixMOutputLayerError.CopyMatrix ( cMatrixOutputLayerError.MergeColumnsToColumnVector () );


		CMatrix	cMatrixJE;
		cMatrixJE.CopyMatrix ( (cMatrixZ.Transpose ()) * cMatrixMOutputLayerError );


		CMatrix cMatrixJJ;
		cMatrixJJ.CopyMatrix ( (cMatrixZ.Transpose ()) * cMatrixZ );
		

		// 定义新的输入层到隐含层的权值
		CMatrix cMatrixNewInputToHideWeight;
		// 定义的新的隐含层的阀值
		CMatrix cMatrixNewHideLayerValve;
		// 定义新的隐含层到输出层的权值
		CMatrix cMatrixNewHideToOutputWeight;
		// 定义新的输出层的阀值
		CMatrix cMatrixNewOutputLayerValve;

		// 定义新的误差矩阵
		CMatrix cMatrixNewOutputLayerError(nOutputLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());

		/////////////////////////////////////////////////////////////////
		// the weight value is adjusted

		while(nStep <= MAX_ADJUST_VALUE)
		{
			CMatrix	cMatrixI (cMatrixZ.GetMatrixColNumber (), cMatrixZ.GetMatrixColNumber ());

			cMatrixI.Eye ();


			CMatrix cMatrixDX;
			cMatrixDX.CopyMatrix ( (((cMatrixJJ + cMatrixI * nStep).Inverse ()) * cMatrixJE) * (-1.0) );
				

			/////////////////////////////////////////////////////////////////////////
			// 拆分cMatrixDX矩阵

			unsigned int nIndex = 0;
				
			// 得到输入层到隐含层的权值的修正量
			CMatrix cMatrixInputToHideWeightChange(nHideLayerNumber, nInputLayerNumber);

			cMatrixDX.CopySubMatrixFromVector (cMatrixInputToHideWeightChange, nIndex);

			nIndex += nHideLayerNumber * nInputLayerNumber;

			// 得到隐含层阀值的修正量
			CMatrix cMatrixHideLayerValveChange (nHideLayerNumber, (unsigned int)1);
				
			cMatrixDX.CopySubMatrixFromVector (cMatrixHideLayerValveChange, nIndex);


			nIndex += nHideLayerNumber;

			// 得到隐含层到输出层的权值的修正量
			CMatrix cMatrixHideToOutputWeightChange (nOutputLayerNumber, nHideLayerNumber);
				
			cMatrixDX.CopySubMatrixFromVector (cMatrixHideToOutputWeightChange, nIndex);

			nIndex += nOutputLayerNumber * nHideLayerNumber;

			// 得到输出层阀值的修正值
			CMatrix cMatrixOutputValveChange (nOutputLayerNumber, (unsigned int)1);

			cMatrixDX.CopySubMatrixFromVector (cMatrixOutputValveChange, nIndex);
	
			cMatrixNewInputToHideWeight.CopyMatrix (matrixInputToHideWeightValue + cMatrixInputToHideWeightChange);

			cMatrixNewHideLayerValve.CopyMatrix (matrixHideLayerValveValue + cMatrixHideLayerValveChange);

			cMatrixNewHideToOutputWeight.CopyMatrix (matrixHideToOutputWeightValue + cMatrixHideToOutputWeightChange);

			cMatrixNewOutputLayerValve.CopyMatrix (matrixOutputLayerValveValue + cMatrixOutputValveChange);
		
			// 前向计算
			::LMForwardCalculate (nInputLayerNumber,
								nHideLayerNumber,
								nOutputLayerNumber,
								//bSimulateDataFlag,
							//	nComboFunc,
								matrixDemoDataInput,
								matrixInputLayerValue,
								cMatrixNewInputToHideWeight,
								cMatrixNewHideLayerValve,
								matrixHideLayerOutput,
								cMatrixNewHideToOutputWeight,
								matrixOutputLayerOutput,
								cMatrixNewOutputLayerValve
								);

				
			cMatrixNewOutputLayerError = cMatrixDemoOutput - matrixOutputLayerOutput;
			nSystemErrorNew = 	cMatrixNewOutputLayerError.GetSystemError ();

			if(nSystemErrorNew < nSystemErrorOld)
			{
				break;
			}
			else
			{
				nStep *= 10;
			}
				
		}// End for while loop

		if ( nStep > MAX_ADJUST_VALUE)
		{
			nLoopTimes--;
			return false;
		}

		nStep	*= 0.1;

		// 赋值
		matrixInputToHideWeightValue = cMatrixNewInputToHideWeight;
		matrixHideLayerValveValue = cMatrixNewHideLayerValve;
		matrixHideToOutputWeightValue = cMatrixNewHideToOutputWeight;
		matrixOutputLayerValveValue = cMatrixNewOutputLayerValve;
		cMatrixOutputLayerError = cMatrixNewOutputLayerError;
		nSystemErrorOld = nSystemErrorNew;

		// 显示数据和程序运行状态
		nSystemError = nSystemErrorOld;
		nTrainTimes = nLoopTimes;


		// 显示系统误差
		//char res_buffer[128] = {0};
		//sprintf(res_buffer,"%.16lf", nSystemError);
		//HWND	hwnd = ::GetDlgItem (hWnd, ID_SYSTEM_ERROR);
		//::SetWindowText (hwnd, res_buffer);

	/*	CString	strSystemError;
		strSystemError.Format ("%.16lf", nSystemError);
		LPCTSTR	lpstrSystemError = (LPCTSTR)strSystemError;
		HWND	hwnd = ::GetDlgItem (hWnd, ID_SYSTEM_ERROR);
		::SetWindowText (hwnd, lpstrSystemError);*/
		
	
		// 显示训练次数
	/*	CString strTrainTimes;
		strTrainTimes.Format ("%u", nTrainTimes + 1);
		LPCTSTR lpstrTrainTimes = (LPCTSTR)strTrainTimes;
		hwnd = ::GetDlgItem (hWnd, ID_TRAIN_TIMES);
		::SetWindowText (hwnd, lpstrTrainTimes);*/

	}// End the "for" loop
	
	return true;
}
void WINAPI BPForwardCalculate ( int nInputLayerNumber,
												int nHideLayerNumber,
												int nOutputLayerNumber,
												//bool bSimulateDataFlag,
												//int nComboFunc,
												CMatrix &matrixDemoDataInput,
												CMatrix &matrixInputLayerValue,
												CMatrix &matrixInputToHideWeightValue,
												CMatrix &matrixHideLayerValveValue,
												CMatrix &matrixHideLayerOutput,
												CMatrix &matrixHideToOutputWeightValue,
												CMatrix &matrixOutputLayerOutput,
												CMatrix &matrixOutputLayerValveValue,
												CMatrix &cMatrixExHideLayerValveValue,
												CMatrix &cMatrixExOutputLayerValveValue
											   )
{

//	if(bSimulateDataFlag)
	//{

		CMatrix cMatrixInputLayerValue(matrixDemoDataInput.GetMatrixRowNumber (), nInputLayerNumber);

		// 得到样本的输入值
		matrixDemoDataInput.CopySubMatrix (cMatrixInputLayerValue, (unsigned int)0, (unsigned int)0);

		CMatrix cMatrixTInputLayerValue = cMatrixInputLayerValue.Transpose ();

		matrixInputLayerValue.CopyMatrix (cMatrixTInputLayerValue);


	//}

	/////////////////////////////////////////////////////////////////////////
	// 得到所有样本的隐含层的净输入
	// 构造规则:
	//		1. 样本的数目做为矩阵行数;
	//		2. 单个样本的隐含层的数目做为矩阵的列数;
	//		3. 矩阵元素中的值即为对应的样本的隐含层的净输入:
	//		   由 
	//				cMatrixInputLayerValue * cMatrixInputToHideWeightValue 
	//				+ cMatrixHideLayerValveValue
	//			得到.
	//
		
	//CMatrix cMatrixExHideLayerValveValue;
	//cMatrixExHideLayerValveValue.nncpyi (matrixHideLayerValveValue, matrixDemoDataInput.GetMatrixRowNumber ());

	CMatrix cMatrixHideLayerPureInput(nHideLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());

	cMatrixHideLayerPureInput = matrixInputToHideWeightValue * matrixInputLayerValue;

	cMatrixHideLayerPureInput += cMatrixExHideLayerValveValue;

	/////////////////////////////////////////////////////////////////////
	//	算出所有样本的隐含层的输出
	//	构造规则:
	//		1. 隐含层的输出y与隐含层的输入x的关系可用函数表示
	//			y = f(x)
	//		2. 矩阵的维数和隐含层的净输入矩阵的维数相等
	//

	CMatrix cMatrixHideLayerOutput(nHideLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());


		cMatrixHideLayerOutput = cMatrixHideLayerPureInput.Tansig();
	
	
	matrixHideLayerOutput.CopyMatrix(cMatrixHideLayerOutput);

	/////////////////////////////////////////////////////////////////////
	//	得到所有样本输出层的净输入
	//	构造规则;
	//		1. 样本的数目做为矩阵的行数;
	//		2. 单个样本的输出层的数目做为矩阵的列数;
	//		3. 矩阵中元素的值即为对应样本的输出层的净输入:
	//			由
	//				cMatrixHideLayerOutput * cMatrixHideToOutputWeightValue
	//				+ cMatrixOutputLayerValveValue
	//			得到
	//

	//CMatrix cMatrixExOutputLayerValveValue;
	//cMatrixExOutputLayerValveValue.nncpyi (matrixOutputLayerValveValue, matrixDemoDataInput.GetMatrixRowNumber ());

	CMatrix cMatrixOutputLayerPureInput(nOutputLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());

	cMatrixOutputLayerPureInput = matrixHideToOutputWeightValue * cMatrixHideLayerOutput;
	cMatrixOutputLayerPureInput  += cMatrixExOutputLayerValveValue;


	/////////////////////////////////////////////////////////////////////
	//	算出所有样本的输出层的输出
	//	构造规则:
	//		1. 矩阵的维数与得到的所有样本的输出层的净输入组成的矩阵一样;
	//		2. 输出层的输出y和输出层的输入可用关系式
	//			y = f(x)
	//		表示
	//

	CMatrix cMatrixOutputLayerOutput(nOutputLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());


	matrixOutputLayerOutput.CopyMatrix(cMatrixOutputLayerOutput);

}

/////////////////////////////////////////////////////////////////////////////
//	Back propagation		---->		前向计算(Only for Simulating)	   //
/////////////////////////////////////////////////////////////////////////////
 void WINAPI BPForwardCalculate2( int nInputLayerNumber,
												int nHideLayerNumber,
												int nOutputLayerNumber,
												//bool bSimulateDataFlag,
												//int nComboFunc,
												CMatrix &matrixDemoDataInput,
												CMatrix &matrixInputLayerValue,
												CMatrix &matrixInputToHideWeightValue,
												CMatrix &matrixHideLayerValveValue,
												CMatrix &matrixHideLayerOutput,
												CMatrix &matrixHideToOutputWeightValue,
												CMatrix &matrixOutputLayerOutput,
												CMatrix &matrixOutputLayerValveValue
											   )
{

	//if(bSimulateDataFlag)
	//{

		CMatrix cMatrixInputLayerValue(matrixDemoDataInput.GetMatrixRowNumber (), nInputLayerNumber);

		// 得到样本的输入值
		matrixDemoDataInput.CopySubMatrix (cMatrixInputLayerValue, (unsigned int)0, (unsigned int)0);

		CMatrix cMatrixTInputLayerValue = cMatrixInputLayerValue.Transpose ();

		matrixInputLayerValue.CopyMatrix (cMatrixTInputLayerValue);


	//}

	/////////////////////////////////////////////////////////////////////////
	// 得到所有样本的隐含层的净输入
	// 构造规则:
	//		1. 样本的数目做为矩阵行数;
	//		2. 单个样本的隐含层的数目做为矩阵的列数;
	//		3. 矩阵元素中的值即为对应的样本的隐含层的净输入:
	//		   由 
	//				cMatrixInputLayerValue * cMatrixInputToHideWeightValue 
	//				+ cMatrixHideLayerValveValue
	//			得到.
	//
		
	CMatrix cMatrixExHideLayerValveValue;
	cMatrixExHideLayerValveValue.nncpyi (matrixHideLayerValveValue, matrixDemoDataInput.GetMatrixRowNumber ());

	CMatrix cMatrixHideLayerPureInput(nHideLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());

	cMatrixHideLayerPureInput = matrixInputToHideWeightValue * matrixInputLayerValue;

	cMatrixHideLayerPureInput += cMatrixExHideLayerValveValue;


	/////////////////////////////////////////////////////////////////////
	//	算出所有样本的隐含层的输出
	//	构造规则:
	//		1. 隐含层的输出y与隐含层的输入x的关系可用函数表示
	//			y = f(x)
	//		2. 矩阵的维数和隐含层的净输入矩阵的维数相等
	//

	CMatrix cMatrixHideLayerOutput(nHideLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());

		cMatrixHideLayerOutput = cMatrixHideLayerPureInput.Tansig();
	
	
	matrixHideLayerOutput.CopyMatrix(cMatrixHideLayerOutput);

	/////////////////////////////////////////////////////////////////////
	//	得到所有样本输出层的净输入
	//	构造规则;
	//		1. 样本的数目做为矩阵的行数;
	//		2. 单个样本的输出层的数目做为矩阵的列数;
	//		3. 矩阵中元素的值即为对应样本的输出层的净输入:
	//			由
	//				cMatrixHideLayerOutput * cMatrixHideToOutputWeightValue
	//				+ cMatrixOutputLayerValveValue
	//			得到
	//

	CMatrix cMatrixExOutputLayerValveValue;
	cMatrixExOutputLayerValveValue.nncpyi (matrixOutputLayerValveValue, matrixDemoDataInput.GetMatrixRowNumber ());

	CMatrix cMatrixOutputLayerPureInput(nOutputLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());

	cMatrixOutputLayerPureInput = matrixHideToOutputWeightValue * cMatrixHideLayerOutput;
	cMatrixOutputLayerPureInput  += cMatrixExOutputLayerValveValue;


	/////////////////////////////////////////////////////////////////////
	//	算出所有样本的输出层的输出
	//	构造规则:
	//		1. 矩阵的维数与得到的所有样本的输出层的净输入组成的矩阵一样;
	//		2. 输出层的输出y和输出层的输入可用关系式
	//			y = f(x)
	//		表示
	//

	CMatrix cMatrixOutputLayerOutput(nOutputLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());



	matrixOutputLayerOutput.CopyMatrix(cMatrixOutputLayerOutput);

}

bool  WINAPI BPDemoDataTrainRepeat (	int nInputLayerNumber,
													int nHideLayerNumber,
													int nOutputLayerNumber,
													//bool bSimulateDataFlag,
													//int nComboFunc,
													double nSystemErrorOld,
													double nSystemErrorNew,
													double nSystemErrorLevel,
													double nSystemError,
													double nStep,
													UINT nMaxTrainTimes,
													UINT nTrainTimes,
													//DWORD ID_SYSTEM_ERROR,
													//DWORD ID_TRAIN_TIMES,
													//HWND  hWnd,
													CMatrix &matrixDemoDataInput,
													CMatrix &matrixInputLayerValue,
													CMatrix &matrixInputToHideWeightValue,
													CMatrix &matrixHideLayerValveValue,
													CMatrix &matrixHideLayerOutput,
													CMatrix &matrixHideToOutputWeightValue,
													CMatrix &matrixOutputLayerOutput,
													CMatrix &matrixOutputLayerValveValue
												 )
{
	// 根据BP算法修正nStep的初始值
	nStep = 0.1;

	// 前向计算
	::LMForwardCalculate (nInputLayerNumber,
						nHideLayerNumber,
						nOutputLayerNumber,
						//bSimulateDataFlag,
						//nComboFunc,
						matrixDemoDataInput,
						matrixInputLayerValue,
						matrixInputToHideWeightValue,
						matrixHideLayerValveValue,
						matrixHideLayerOutput,
						matrixHideToOutputWeightValue,
						matrixOutputLayerOutput,
						matrixOutputLayerValveValue
						);
		

	/////////////////////////////////////////////////////////////////////
	//	算出所有样本的输出层的delta矩阵
	//	构造规则:
	//		1. 样本的数目为矩阵的行数;
	//		2. 样本输出层的数目为矩阵的列数;
	//		3. 矩阵中的元素的值y为:
	//				y = (前向计算出的输出层的值 - 样本的输出层的值) .* f'(net)
	//
	CMatrix cMatrixTDemoOutput(matrixDemoDataInput.GetMatrixRowNumber (), nOutputLayerNumber);
		
	// 得到样本中输出层的数据
	matrixDemoDataInput.CopySubMatrix (cMatrixTDemoOutput, (unsigned int)nInputLayerNumber, (unsigned int)0);

	CMatrix cMatrixDemoOutput = cMatrixTDemoOutput.Transpose ();

	// 得到样本中输出层的误差
	CMatrix cMatrixOutputLayerError(nOutputLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());
	cMatrixOutputLayerError = cMatrixDemoOutput - matrixOutputLayerOutput;

	nSystemErrorOld = cMatrixOutputLayerError.GetSystemError ();

	for(int nLoopTimes=1; nLoopTimes < nMaxTrainTimes; nLoopTimes++)	
	{
		if(nSystemErrorOld < nSystemErrorLevel)
		{
			nLoopTimes--;
			break;
		}

		// 求输出层的delta值
		// 注意: 此处'/' 是 '点乘'!!!
		CMatrix	cMatrixOutputLayerDelta (nOutputLayerNumber, matrixDemoDataInput.GetMatrixRowNumber());
		cMatrixOutputLayerDelta = (matrixOutputLayerOutput - matrixOutputLayerOutput / matrixOutputLayerOutput) / cMatrixOutputLayerError;
	
		CMatrix cMatrixTHideToOutputWeightValue (matrixHideToOutputWeightValue.GetMatrixColNumber(), matrixHideToOutputWeightValue.GetMatrixRowNumber());
		cMatrixTHideToOutputWeightValue = matrixHideToOutputWeightValue.Transpose();

		// 求隐含层的delta值
		// 注意: 此处'/' 是 '点乘'!!!
		CMatrix cMatrixHideLayerDelta;
		cMatrixHideLayerDelta.CopyMatrix ( (matrixHideLayerOutput - (matrixHideLayerOutput / matrixHideLayerOutput)) / ( cMatrixTHideToOutputWeightValue * cMatrixOutputLayerDelta) );
				
		// 定义新的输入层到隐含层的权值
		CMatrix cMatrixNewInputToHideWeight (matrixInputToHideWeightValue.GetMatrixRowNumber (), matrixInputToHideWeightValue.GetMatrixColNumber ());
		// 定义的新的隐含层的阀值
		CMatrix cMatrixNewHideLayerValve (nHideLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());
		// 定义新的隐含层到输出层的权值
		CMatrix cMatrixNewHideToOutputWeight (matrixHideToOutputWeightValue.GetMatrixRowNumber (), matrixHideToOutputWeightValue.GetMatrixColNumber ());
		// 定义新的输出层的阀值
		CMatrix cMatrixNewOutputLayerValve (nOutputLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());
		// 定义新的误差矩阵
		CMatrix cMatrixNewOutputLayerError(nOutputLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());


		// 权值和阀值调整
		cMatrixNewHideToOutputWeight = cMatrixOutputLayerDelta * (matrixHideLayerOutput.Transpose ()) * (nStep);
		cMatrixNewOutputLayerValve = cMatrixOutputLayerDelta;

		cMatrixNewInputToHideWeight = cMatrixHideLayerDelta * (matrixInputLayerValue.Transpose ()) * (nStep);
		cMatrixNewHideLayerValve = cMatrixHideLayerDelta;

		// 赋值
		matrixInputToHideWeightValue += cMatrixNewInputToHideWeight;

		CMatrix cMatrixExHideLayerValveValue;
		cMatrixExHideLayerValveValue.nncpyi (matrixHideLayerValveValue, matrixDemoDataInput.GetMatrixRowNumber ());
		cMatrixExHideLayerValveValue += cMatrixNewHideLayerValve;

		matrixHideToOutputWeightValue += cMatrixNewHideToOutputWeight;

		CMatrix cMatrixExOutputLayerValveValue;
		cMatrixExOutputLayerValveValue.nncpyi (matrixOutputLayerValveValue, matrixDemoDataInput.GetMatrixRowNumber ());
		cMatrixExOutputLayerValveValue += cMatrixNewOutputLayerValve;

		// 前向计算
		::BPForwardCalculate (nInputLayerNumber,
							nHideLayerNumber,
							nOutputLayerNumber,
							//bSimulateDataFlag,
							//nComboFunc,
							matrixDemoDataInput,
							matrixInputLayerValue,
							matrixInputToHideWeightValue,
							matrixHideLayerValveValue,
							matrixHideLayerOutput,
							matrixHideToOutputWeightValue,
							matrixOutputLayerOutput,
							matrixOutputLayerValveValue,
							cMatrixExHideLayerValveValue,
							cMatrixExOutputLayerValveValue
							);


		cMatrixNewOutputLayerError = cMatrixDemoOutput - matrixOutputLayerOutput;
		nSystemErrorNew = 	cMatrixNewOutputLayerError.GetSystemError ();

		cMatrixOutputLayerError = cMatrixNewOutputLayerError;

		if(nSystemErrorNew < nSystemErrorOld)
		{
			nSystemErrorOld = nSystemErrorNew;
		}
		else
		{
			nStep *= -0.1;
		}

		// 显示数据和程序运行状态
		nSystemError = nSystemErrorOld;
		nTrainTimes = nLoopTimes;


		// 显示系统误差
	/*	CString	strSystemError;
		strSystemError.Format ("%lf", nSystemError);
		LPCTSTR	lpstrSystemError = (LPCTSTR)strSystemError;
		HWND	hwnd = ::GetDlgItem (hWnd, ID_SYSTEM_ERROR);
		::SetWindowText (hwnd, lpstrSystemError);*/
	
		// 显示训练次数
		/*CString strTrainTimes;
		strTrainTimes.Format ("%u", nTrainTimes + 1);
		LPCTSTR lpstrTrainTimes = (LPCTSTR)strTrainTimes;
		hwnd = ::GetDlgItem (hWnd, ID_TRAIN_TIMES);
		::SetWindowText (hwnd, lpstrTrainTimes);*/

	}// End the "for" loop
	
	return true;
}



