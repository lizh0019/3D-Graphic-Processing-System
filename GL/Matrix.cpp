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
	// ��̬�����ά����
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

	// �Զ��������ֵ
	m_nRow	= nRow;
	m_nCol	= nCol;
	m_pTMatrix = tMatrix;

}

CMatrix::CMatrixone(unsigned int nRow,unsigned int nCol)
{
	// ��̬�����ά����
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

	// �Զ��������ֵ
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
	// Ҫ���������ӵ�����: ������Ŀ���!
	if(m_nRow != cMatrixB.m_nRow || m_nCol != cMatrixB.m_nCol )
	{
		::AfxMessageBox (TEXT("ִ����ӵ���������ά�������!"),MB_OK | MB_ICONERROR);
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
	// Ҫ���������ӵ�����: ��������Ŀ���!
	if(m_nRow != cMatrixB.m_nRow || m_nCol != cMatrixB.m_nCol )
	{
		::AfxMessageBox (TEXT("ִ���������������ά�������!"),MB_OK | MB_ICONERROR);
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
		::AfxMessageBox (TEXT("ִ����˵���������ά����������˵�����!"),MB_OK | MB_ICONERROR);
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
		::AfxMessageBox(TEXT("�Ⱥ��������ߵľ����ά�������!"),MB_OK | MB_ICONERROR);
		return *this;	// return invalid value
	}	

	// ��������ֵ
	m_nRow = cMatrixB.m_nRow ;
	m_nCol = cMatrixB.m_nCol ;
	m_pTMatrix = cMatrixB.m_pTMatrix ;

	// ��ֵ����
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
	
	// ��������ֵ
	m_nRow = cMatrixB.m_nRow ;
	m_nCol = cMatrixB.m_nCol ;
	m_pTMatrix = cMatrixB.m_pTMatrix ;

	// ��ֵ����
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
		//printf("����!ִ����ӵ���������ά�������!\n");
		::AfxMessageBox (TEXT("����������߾����ά�������!"),MB_OK | MB_ICONERROR);
		return *this;	// return invalid value
	}
	
	// ��ֵ����
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

	// �Ծ���ֵ
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
	// �ο���Ŀ: �������ֵ���� --->ʩ���� �¹�֦
	/////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////
	// �ж��Ƿ��ǿ�����:
	//		������һ���Ƿ���!!!

	if ( m_nRow != m_nCol)
	{
		//printf("����!����������������,�Ƿǿ�����!\n");
		::AfxMessageBox (TEXT("����������������,�Ƿǿ�����!"),MB_OK | MB_ICONERROR);
	}

	CMatrix cMatrix = *this;

	//***********************************************************************
	// ˼·:(�ǳ���˼ά!)
	//		��̬������������(2*m_nCol)���洢ÿ�ν������������ֵ
	//		������û���н�������¼��������,
	//		1.û�����н����������������,��SwapMatrixRow()������
	//		��⵽����ֵ�����������,��SwapMatrixCol()������Ҳһ��,
	//		��⵽����ֵ�����������,��ռ��ϵͳ��Դ;
	//		2.����ȵľͽ���
	//***********************************************************************
    
	//	�����ڴ�
	int *pIntArray = new int [2*m_nCol];

	// nSetp --- Լ������,����չ��	
	for(unsigned int k=0; k < cMatrix.m_nCol; k++)
	{
		/////////////////////////////////////////////////////////////////////
		// �����н��� ---> ��Ϸ����:
		// Ϊ��֤������̵���ֵ�ȶ���,�ڵ�k��Լ��ʱ,����{a(ik)}|i=k->n��ѡ��
		// ģ�������ΪԼ����Ԫ��,������������Ӧ����

		// �����Ԫ��
		double nMaxElement = cMatrix.m_pTMatrix [k][k];
		// �����Ԫ�����ڵ�����
		unsigned int nMainRow = k;

		for(unsigned int nCount = k+1; nCount < cMatrix.m_nCol; nCount++)
		{
			if( fabs(nMaxElement) < fabs(cMatrix.m_pTMatrix [nCount][k]) )
			{
				nMaxElement = cMatrix.m_pTMatrix [nCount][k];
				nMainRow = nCount;
			}
		}

		// ������������������������
		pIntArray [2*k] = k;
		pIntArray [2*k+1] = nMainRow;
	

		// ������
		cMatrix.SwapMatrixRow(k,nMainRow);

		//Display();

		//	�ж��Ƿ��ǿ�����
		if(cMatrix.m_pTMatrix [k][k] == 0)
		{
			//printf("����!�˾���Ϊ�ǿ�����!\n");
			::AfxMessageBox (TEXT("�˾���Ϊ�ǿ�����,û�������!"),MB_OK | MB_ICONERROR);
		}

		cMatrix.m_pTMatrix [k][k] = 1/(cMatrix.m_pTMatrix [k][k]);
		

		// ������
		for(unsigned int i=0; i < cMatrix.m_nRow; i++)
		{	
			if( i != k)
				cMatrix.m_pTMatrix [i][k] = -(cMatrix.m_pTMatrix [k][k]) * (cMatrix.m_pTMatrix [i][k]); 
			
			//int nTempValue = m_pTMatrix [i][k];
			
		}
		
		//printf("\n");

		// Լ��������
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

		// ������
		for(unsigned int j=0; j < cMatrix.m_nCol; j++)
		{
			if( j != k)
				cMatrix.m_pTMatrix [k][j] = (cMatrix.m_pTMatrix [k][k]) * (cMatrix.m_pTMatrix [k][j]);

		}

	}

	
	/////////////////////////////////////////////////////////////////////////
	// �����н��� ---> �Խ����к�ľ�������н��� ---> ��ԭ����
	// ��Ϸ����:
	// ����ʼ�����н��е��н��� ---> �������Ӧ���н������л�ԭ,���ɵõ������
	// �����

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
		::AfxMessageBox (TEXT("���ܴ�Ҫ��ȡ���ݵ��ļ�!"),MB_OK | MB_ICONERROR);
		dataFile.Close ();
		return FALSE;
	}

	// �����洢��ȡ�ı��ļ���һ�е�����
	CString	strData;				 
	// ������¼�ı��ļ���һ���ж���������?
	unsigned int	nRow = 0;		

	/////////////////////////////////////////////////////////////////////////
	// Step 1: �õ��ļ���������Ŀ�������ı��ļ���������Ŀ�����ö���(����)����
	// ����Ŀ
	//

	while(dataFile.ReadString (strData) != FALSE)
	{
		++nRow;
	}

	// �����ı��ļ������ݵ��������ö���(����)������
	m_nRow = nRow;
	SetMatrixRowNumber(m_nRow);

	// ���¶�λ��ǰ�ļ�ָ�뵽�ļ���ͷ
	dataFile.SeekToBegin ();

	dataFile.ReadString (strData);
	strData.TrimLeft ();
	strData.TrimRight ();

	TCHAR	SPACE_CHARACTER = ' ';
	// ������¼�ı��ļ���һ���ж�����?
	unsigned int	nCol = 0;						
	
	// �ո�����ַ����е�λ������
	int nIndex = 0;

	do
	{
		nIndex = strData.Find (SPACE_CHARACTER);

		++nCol;

		// ��ȡ�ַ��������ַ���,����ȡһ��double��ʵ������
		//CString strDoubleNumber = strData.Left (nIndex);	

		// ���ַ���ת��Ϊdouble��ʵ��
		//double RealNumber = atof(strDoubleNumber);

		//int nTempNumber = strData.GetLength ();

		strData = strData.Right (strData.GetLength () - nIndex -1);

		// ȥ������Ŀո�
		strData.TrimLeft ();

		// Use for debugging
		//int nTempNum = strData.GetLength ();

	}while(nIndex != -1);

	// �����ı��ļ������ݵ��������ö���(����)������
	m_nCol = nCol;
	SetMatrixColNumber(m_nCol);

	// End of Getting the Rows and Cols of the Text File
	/////////////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////////////////////////////////
	// Step 2: �����ı��ļ��е����ݶԾ���ֵ,�����ÿ�е������Ƿ�͵�һ�е�
	// �������,�������ʾ������Ϣ
	//

	// ���¶�λ��ǰ�ļ�ָ�뵽�ļ���ͷ
	dataFile.SeekToBegin ();

	// �Ծ����е�Ԫ��װ���ı��ļ�������
	for(unsigned int i=0; i < m_nRow; i++)
	{
		dataFile.ReadString (strData);
		strData.TrimLeft ();
		strData.TrimRight ();

		// ��֤ÿ�е������Ƿ����һ�е��������
		unsigned int nVerifyColNum = 0;

		do
		{
			nIndex = strData.Find (SPACE_CHARACTER);

			++nVerifyColNum;

			if(nIndex != -1)
			{
				// ��ȡ�ַ��������ַ���,����ȡһ��double��ʵ������
				CString strDoubleNumber = strData.Left (nIndex);	

				// ���ַ���ת��Ϊdouble��ʵ��
				double RealNumber = atof(strDoubleNumber);

				m_pTMatrix [i][nVerifyColNum - 1] = RealNumber;
				//int nTempNumber = strData.GetLength ();

				strData = strData.Right (strData.GetLength () - nIndex -1);

				// ȥ������Ŀո�
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

			CString strMessage =  CString(TEXT("�ı��ļ���")) + strRowNumber + CString(TEXT("��һ����")) + strVerifyColNumber + CString(TEXT("��,���һ���е�����")) + strColNumber + CString(TEXT("�����!"));

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
		::AfxMessageBox(TEXT("���ܴ����ļ�!"),MB_OK | MB_ICONERROR);
		dataFile.Close ();
		return false;
	}
	
	dataFile.SeekToEnd ();

	// ������(����)�е�����д���ļ�
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
		::AfxMessageBox (TEXT("�˾���������������!����ת��Ϊ��λ��!"),MB_OK | MB_ICONERROR);
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
		::AfxMessageBox (TEXT("�����ľ�����������!"),MB_OK | MB_ICONERROR);
		return;
	}

	if((m_nRow - nIndex) < (cMatrix.m_nRow * cMatrix.m_nCol))
	{
		::AfxMessageBox (TEXT("��������Ŀռ���������!"),MB_OK | MB_ICONERROR);
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
		::AfxMessageBox (TEXT("�����������������,����������!"),MB_OK | MB_ICONERROR);
		return;
	}

	// Verify whether the number of the object element is enough to be copyed

	if((m_nRow - nIndex) < (cMatrix.m_nRow * cMatrix.m_nCol))
	{
		::AfxMessageBox (TEXT("�����е�Ԫ����������!"),MB_OK | MB_ICONERROR);
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
		::AfxMessageBox (TEXT("�����������������,����������!"),MB_OK | MB_ICONERROR);
		return;
	}

	// Verify whether the number of the object element is enough to be copyed
	if((m_nRow - nIndex) < cMatrix.m_nCol )
	{
		::AfxMessageBox (TEXT("�����Ԫ����������!"),MB_OK | MB_ICONERROR);
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
		::AfxMessageBox (TEXT("�����ľ�����������!"),MB_OK | MB_ICONERROR);
		return;
	}

	if((m_nRow - nIndex) < cMatrix.m_nCol)
	{
		::AfxMessageBox (TEXT("��������Ŀռ���������!"),MB_OK | MB_ICONERROR);
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
// �Ծ��������е�Ԫ�ؽ���һ�η����Ա任:
//		�任���ֵy��任ǰ��ֵ�Ĺ�ϵ��:
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
// �Ծ��������е�Ԫ�ؽ���һ�η����Ա任:
//		�任���ֵy��任ǰ��ֵ�Ĺ�ϵ��:
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
// �Ծ��������е�Ԫ�ؽ���һ�η����Ա任:
//		�任���ֵy��任ǰ��ֵ�Ĺ�ϵ��:
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
// �Ծ��������е�Ԫ�ؽ���һ�η����Ա任:
//		�任���ֵy��任ǰ��ֵ�Ĺ�ϵ��:
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
// �Ծ��������е�Ԫ�ؽ���һ�η����Ա任:
//		�任���ֵy��任ǰ��ֵ�Ĺ�ϵ��:
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
// �Ծ��������е�Ԫ�ؽ���һ�η����Ա任:
//		�任���ֵy��任ǰ��ֵ�Ĺ�ϵ��:
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

	// �����ڴ�
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

	// ���ݿռ�����ڴ�
	m_pTMatrix.resize (m_nRow);
	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			m_pTMatrix[i].resize (m_nCol);
			m_pTMatrix[i][j] = (double) 0;
		}
	}

	// �Ծ���ֵ
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
		::AfxMessageBox (TEXT("�������ľ���ά��С��Ҫ�����ľ�������Ҫ��ά��!"),MB_OK | MB_ICONERROR);
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
		::AfxMessageBox (TEXT("�������ľ�����������!!!"),MB_OK | MB_ICONERROR);
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

	// ���ݿռ�����ڴ�
	m_pTMatrix.resize (m_nRow);
	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			m_pTMatrix[i].resize (m_nCol);
			m_pTMatrix[i][j] = (BYTE) 0;
		}
	}

	// ��ֵ
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

	// ���ݿռ�����ڴ�
	m_pTMatrix.resize (m_nRow);
	for(unsigned int i=0; i < m_nRow; i++)
	{
		for(unsigned int j=0; j < m_nCol; j++)
		{
			m_pTMatrix[i].resize (m_nCol);
			m_pTMatrix[i][j] = (BYTE) 0;
		}
	}

	// ������ֵ
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
		::AfxMessageBox (TEXT("���������ά�������,����������˵�����!"),MB_OK | MB_ICONERROR);
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
	//	�������
	if( cMatrixA.GetMatrixRowNumber () != cMatrixB.GetMatrixRowNumber () )
	{
		::AfxMessageBox (TEXT("����ϲ���������������������!"),MB_OK | MB_ICONERROR);
	
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
	// ���������Ԫ�صľ���
	// �������:
	//		1. ������Ŀ��Ϊ���������;
	//		2. ������������������Ŀ��Ϊ���������;
	//		3. �����е�Ԫ�ؼ�Ϊ��Ӧ��������ֵ
	//

	CMatrix cMatrixInputLayerValue(matrixDemoDataInput.GetMatrixRowNumber (), nInputLayerNumber);

	// �õ�����������ֵ
	matrixDemoDataInput.CopySubMatrix (cMatrixInputLayerValue,(unsigned int)0,(unsigned int)0);

	CMatrix cMatrixTInputLayerValue = cMatrixInputLayerValue.Transpose ();

	matrixInputLayerValue.CopyMatrix (cMatrixTInputLayerValue);


	/////////////////////////////////////////////////////////////////////////
	// ����Ȩֵ���� -----> �ɵ��������������������֮���Ȩֵ��ΪԪ�����
	// �������:
	//		1. ������������������Ŀ��Ϊ��������;
	//		2. �������������������Ŀ��Ϊ���������;
	//		3. �Ծ����е�Ԫ�ؽ��������ʼ��,ֵ��(-1,1)֮��;
	//		4. �������������������������Ŀ����ȵ�;
	//		5. ��������ʹ�õ���ͬһ��Ȩֵ����.
	//

	CMatrix	cMatrixInputToHideWeightValue(nHideLayerNumber, nInputLayerNumber);

	// �����ʼ��������Ԫ�ص�ֵ
	cMatrixInputToHideWeightValue.RandomInitialize ();

	matrixInputToHideWeightValue.CopyMatrix (cMatrixInputToHideWeightValue);


	/////////////////////////////////////////////////////////////////////
	//	��������������ķ�ֵ����
	//	�������:
	//		1. ��������Ŀ��Ϊ��������;
	//		2. �������������������Ŀ��Ϊ���������;
	//		3. �����ʼ�������е�Ԫ��,ʹֵ��(-1,1)֮��;
	//		4. ������ÿ�е����ݶ��͵�һ���������Ӧ��λ�����.
	//

	CMatrix cMatrixHideLayerValveValue(nHideLayerNumber,(unsigned int)1);

	// �����ʼ��������Ԫ�ص�ֵ
	cMatrixHideLayerValveValue.RandomInitialize ();

	matrixHideLayerValveValue.CopyMatrix(cMatrixHideLayerValveValue);


	/////////////////////////////////////////////////////////////////////
	//	����Ȩֵ����	----->	�ɵ����������������������֮��Ȩֵ��ΪԪ��
	//							���
	//	�������:
	//		1. �������������������Ŀ��Ϊ���������;
	//		2. ������������������Ŀ��Ϊ���������;
	//		3. �Ծ����е�Ԫ�ؽ��������ʼ��,ֵ��(-1,1)֮��;
	//		4. �������������������������Ŀ����ȵ�;
	//		5. ��������ʹ�õ���ͬһ��Ȩֵ����.
	//

	CMatrix cMatrixHideToOutputWeightValue(nOutputLayerNumber, nHideLayerNumber);

	// �Ծ����Ԫ�������ʼ��
	cMatrixHideToOutputWeightValue.RandomInitialize ();

	matrixHideToOutputWeightValue.CopyMatrix (cMatrixHideToOutputWeightValue);


	/////////////////////////////////////////////////////////////////////
	//	���������������ķ�ֵ����
	//	�������:
	//		1. ��������Ŀ��Ϊ���������;
	//		2. ������������������Ŀ��Ϊ���������;
	//		3. �����ʼ�������е�Ԫ��,ʹֵ��(-1,1)֮��;
	//		4. ������ÿ�е����ݶ��͵�һ���������Ӧ��λ�����.
	//

	CMatrix cMatrixOutputLayerValveValue(nOutputLayerNumber,(unsigned int)1);

	// �����ʼ��������Ԫ�ص�ֵ
	cMatrixOutputLayerValveValue.RandomInitialize ();

	matrixOutputLayerValveValue.CopyMatrix(cMatrixOutputLayerValveValue);


}
/////////////////////////////////////////////////////////////////////////////
//	Levenberg-Marquart		---->		ǰ�����						   //
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

		// �õ�����������ֵ
		matrixDemoDataInput.CopySubMatrix (cMatrixInputLayerValue, (unsigned int)0, (unsigned int)0);

		CMatrix cMatrixTInputLayerValue = cMatrixInputLayerValue.Transpose ();

//cMatrixTInputLayerValue =cMatrixTInputLayerValue.Tansig();
		matrixInputLayerValue.CopyMatrix (cMatrixTInputLayerValue);





	/////////////////////////////////////////////////////////////////////////
	// �õ�����������������ľ�����
	// �������:
	//		1. ��������Ŀ��Ϊ��������;
	//		2. �������������������Ŀ��Ϊ���������;
	//		3. ����Ԫ���е�ֵ��Ϊ��Ӧ��������������ľ�����:
	//		   �� 
	//				cMatrixInputLayerValue * cMatrixInputToHideWeightValue 
	//				+ cMatrixHideLayerValveValue
	//			�õ�.
	//
		
	CMatrix cMatrixExHideLayerValveValue;
	cMatrixExHideLayerValveValue.nncpyi (matrixHideLayerValveValue, matrixDemoDataInput.GetMatrixRowNumber ());

	CMatrix cMatrixHideLayerPureInput(nHideLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());

	cMatrixHideLayerPureInput = matrixInputToHideWeightValue * matrixInputLayerValue;

	cMatrixHideLayerPureInput += cMatrixExHideLayerValveValue;


	/////////////////////////////////////////////////////////////////////
	//	�����������������������
	//	�������:
	//		1. ����������y�������������x�Ĺ�ϵ���ú�����ʾ
	//			y = f(x)
	//		2. �����ά����������ľ���������ά�����
	//

	CMatrix cMatrixHideLayerOutput(nHideLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());

	
		cMatrixHideLayerOutput = cMatrixHideLayerPureInput.Tansig();
	

	
	matrixHideLayerOutput.CopyMatrix(cMatrixHideLayerOutput);

	/////////////////////////////////////////////////////////////////////
	//	�õ��������������ľ�����
	//	�������;
	//		1. ��������Ŀ��Ϊ���������;
	//		2. ������������������Ŀ��Ϊ���������;
	//		3. ������Ԫ�ص�ֵ��Ϊ��Ӧ�����������ľ�����:
	//			��
	//				cMatrixHideLayerOutput * cMatrixHideToOutputWeightValue
	//				+ cMatrixOutputLayerValveValue
	//			�õ�
	//

	CMatrix cMatrixExOutputLayerValveValue;
	cMatrixExOutputLayerValveValue.nncpyi (matrixOutputLayerValveValue, matrixDemoDataInput.GetMatrixRowNumber ());

	CMatrix cMatrixOutputLayerPureInput(nOutputLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());

	cMatrixOutputLayerPureInput = matrixHideToOutputWeightValue * cMatrixHideLayerOutput;
	cMatrixOutputLayerPureInput  += cMatrixExOutputLayerValveValue;


	/////////////////////////////////////////////////////////////////////
	//	����������������������
	//	�������:
	//		1. �����ά����õ������������������ľ�������ɵľ���һ��;
	//		2. ���������y��������������ù�ϵʽ
	//			y = f(x)
	//		��ʾ
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
	//	�������������������delta����
	//	�������:
	//		1. ��������ĿΪ���������;
	//		2. ������������ĿΪ���������;
	//		3. �����е�Ԫ�ص�ֵyΪ:
	//				y = -(ǰ��������������ֵ - ������������ֵ) * f'(net)
	//


	CMatrix cMatrixTDemoOutput(matrixDemoDataInput.GetMatrixRowNumber (), nOutputLayerNumber);
		
	// �õ�����������������????
	matrixDemoDataInput.CopySubMatrix (cMatrixTDemoOutput, (unsigned int)nInputLayerNumber, (unsigned int)0);

	CMatrix cMatrixDemoOutput = cMatrixTDemoOutput.Transpose ();

	// �õ����������������
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
		// ע��: �˴�'/' �� '���'!!!
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
		// ע��: �˴�'/' �� '���'!!!
		cMatrixExHideLayerDelta.CopyMatrix ( (1 - (cMatrixExHideLayerOutput / cMatrixExHideLayerOutput)) / ( cMatrixTHideToOutputWeightValue * cMatrixExOutputLayerDelta) );
			//???	

		CMatrix	cMatrixExInputLayerValue;
		cMatrixExInputLayerValue.nncpyi (matrixInputLayerValue, nOutputLayerNumber);

		
		CMatrix	cMatrixJ11;
		cMatrixJ11.nncpy (cMatrixExHideLayerDelta.Transpose(), cMatrixExInputLayerValue.GetMatrixRowNumber ());

		CMatrix cMatrixJ12;
		cMatrixJ12.nncpyi(cMatrixExInputLayerValue.Transpose (), cMatrixExHideLayerDelta.GetMatrixRowNumber());


		CMatrix	cMatrixJ1;
		// ע��: �˴�'/' �� '���'!!!
		cMatrixJ1.CopyMatrix (cMatrixJ11 / cMatrixJ12);


		CMatrix cMatrixJ21;
		cMatrixJ21.nncpy (cMatrixExOutputLayerDelta.Transpose (), cMatrixExHideLayerOutput.GetMatrixRowNumber ());

		CMatrix cMatrixJ22;
		cMatrixJ22.nncpyi (cMatrixExHideLayerOutput.Transpose (), cMatrixExOutputLayerDelta.GetMatrixRowNumber ());

		CMatrix cMatrixJ2;
		// ע��: �˴�'/' �� '���'!!!
		cMatrixJ2.CopyMatrix (cMatrixJ21 / cMatrixJ22);

		
		CMatrix cMatrixZ;
		cMatrixZ.CopyMatrix (MergeMatrix(MergeMatrix(MergeMatrix(cMatrixJ1, cMatrixExHideLayerDelta.Transpose ()), cMatrixJ2), cMatrixExOutputLayerDelta.Transpose ()));


		CMatrix cMatrixMOutputLayerError;
		cMatrixMOutputLayerError.CopyMatrix ( cMatrixOutputLayerError.MergeColumnsToColumnVector () );


		CMatrix	cMatrixJE;
		cMatrixJE.CopyMatrix ( (cMatrixZ.Transpose ()) * cMatrixMOutputLayerError );


		CMatrix cMatrixJJ;
		cMatrixJJ.CopyMatrix ( (cMatrixZ.Transpose ()) * cMatrixZ );
		

		// �����µ�����㵽�������Ȩֵ
		CMatrix cMatrixNewInputToHideWeight;
		// ������µ�������ķ�ֵ
		CMatrix cMatrixNewHideLayerValve;
		// �����µ������㵽������Ȩֵ
		CMatrix cMatrixNewHideToOutputWeight;
		// �����µ������ķ�ֵ
		CMatrix cMatrixNewOutputLayerValve;

		// �����µ�������
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
			// ���cMatrixDX����

			unsigned int nIndex = 0;
				
			// �õ�����㵽�������Ȩֵ��������
			CMatrix cMatrixInputToHideWeightChange(nHideLayerNumber, nInputLayerNumber);

			cMatrixDX.CopySubMatrixFromVector (cMatrixInputToHideWeightChange, nIndex);

			nIndex += nHideLayerNumber * nInputLayerNumber;

			// �õ������㷧ֵ��������
			CMatrix cMatrixHideLayerValveChange (nHideLayerNumber, (unsigned int)1);
				
			cMatrixDX.CopySubMatrixFromVector (cMatrixHideLayerValveChange, nIndex);


			nIndex += nHideLayerNumber;

			// �õ������㵽������Ȩֵ��������
			CMatrix cMatrixHideToOutputWeightChange (nOutputLayerNumber, nHideLayerNumber);
				
			cMatrixDX.CopySubMatrixFromVector (cMatrixHideToOutputWeightChange, nIndex);

			nIndex += nOutputLayerNumber * nHideLayerNumber;

			// �õ�����㷧ֵ������ֵ
			CMatrix cMatrixOutputValveChange (nOutputLayerNumber, (unsigned int)1);

			cMatrixDX.CopySubMatrixFromVector (cMatrixOutputValveChange, nIndex);
	
			cMatrixNewInputToHideWeight.CopyMatrix (matrixInputToHideWeightValue + cMatrixInputToHideWeightChange);

			cMatrixNewHideLayerValve.CopyMatrix (matrixHideLayerValveValue + cMatrixHideLayerValveChange);

			cMatrixNewHideToOutputWeight.CopyMatrix (matrixHideToOutputWeightValue + cMatrixHideToOutputWeightChange);

			cMatrixNewOutputLayerValve.CopyMatrix (matrixOutputLayerValveValue + cMatrixOutputValveChange);
		
			// ǰ�����
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

		// ��ֵ
		matrixInputToHideWeightValue = cMatrixNewInputToHideWeight;
		matrixHideLayerValveValue = cMatrixNewHideLayerValve;
		matrixHideToOutputWeightValue = cMatrixNewHideToOutputWeight;
		matrixOutputLayerValveValue = cMatrixNewOutputLayerValve;
		cMatrixOutputLayerError = cMatrixNewOutputLayerError;
		nSystemErrorOld = nSystemErrorNew;

		// ��ʾ���ݺͳ�������״̬
		nSystemError = nSystemErrorOld;
		nTrainTimes = nLoopTimes;


		// ��ʾϵͳ���
		//char res_buffer[128] = {0};
		//sprintf(res_buffer,"%.16lf", nSystemError);
		//HWND	hwnd = ::GetDlgItem (hWnd, ID_SYSTEM_ERROR);
		//::SetWindowText (hwnd, res_buffer);

	/*	CString	strSystemError;
		strSystemError.Format ("%.16lf", nSystemError);
		LPCTSTR	lpstrSystemError = (LPCTSTR)strSystemError;
		HWND	hwnd = ::GetDlgItem (hWnd, ID_SYSTEM_ERROR);
		::SetWindowText (hwnd, lpstrSystemError);*/
		
	
		// ��ʾѵ������
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

		// �õ�����������ֵ
		matrixDemoDataInput.CopySubMatrix (cMatrixInputLayerValue, (unsigned int)0, (unsigned int)0);

		CMatrix cMatrixTInputLayerValue = cMatrixInputLayerValue.Transpose ();

		matrixInputLayerValue.CopyMatrix (cMatrixTInputLayerValue);


	//}

	/////////////////////////////////////////////////////////////////////////
	// �õ�����������������ľ�����
	// �������:
	//		1. ��������Ŀ��Ϊ��������;
	//		2. �������������������Ŀ��Ϊ���������;
	//		3. ����Ԫ���е�ֵ��Ϊ��Ӧ��������������ľ�����:
	//		   �� 
	//				cMatrixInputLayerValue * cMatrixInputToHideWeightValue 
	//				+ cMatrixHideLayerValveValue
	//			�õ�.
	//
		
	//CMatrix cMatrixExHideLayerValveValue;
	//cMatrixExHideLayerValveValue.nncpyi (matrixHideLayerValveValue, matrixDemoDataInput.GetMatrixRowNumber ());

	CMatrix cMatrixHideLayerPureInput(nHideLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());

	cMatrixHideLayerPureInput = matrixInputToHideWeightValue * matrixInputLayerValue;

	cMatrixHideLayerPureInput += cMatrixExHideLayerValveValue;

	/////////////////////////////////////////////////////////////////////
	//	�����������������������
	//	�������:
	//		1. ����������y�������������x�Ĺ�ϵ���ú�����ʾ
	//			y = f(x)
	//		2. �����ά����������ľ���������ά�����
	//

	CMatrix cMatrixHideLayerOutput(nHideLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());


		cMatrixHideLayerOutput = cMatrixHideLayerPureInput.Tansig();
	
	
	matrixHideLayerOutput.CopyMatrix(cMatrixHideLayerOutput);

	/////////////////////////////////////////////////////////////////////
	//	�õ��������������ľ�����
	//	�������;
	//		1. ��������Ŀ��Ϊ���������;
	//		2. ������������������Ŀ��Ϊ���������;
	//		3. ������Ԫ�ص�ֵ��Ϊ��Ӧ�����������ľ�����:
	//			��
	//				cMatrixHideLayerOutput * cMatrixHideToOutputWeightValue
	//				+ cMatrixOutputLayerValveValue
	//			�õ�
	//

	//CMatrix cMatrixExOutputLayerValveValue;
	//cMatrixExOutputLayerValveValue.nncpyi (matrixOutputLayerValveValue, matrixDemoDataInput.GetMatrixRowNumber ());

	CMatrix cMatrixOutputLayerPureInput(nOutputLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());

	cMatrixOutputLayerPureInput = matrixHideToOutputWeightValue * cMatrixHideLayerOutput;
	cMatrixOutputLayerPureInput  += cMatrixExOutputLayerValveValue;


	/////////////////////////////////////////////////////////////////////
	//	����������������������
	//	�������:
	//		1. �����ά����õ������������������ľ�������ɵľ���һ��;
	//		2. ���������y��������������ù�ϵʽ
	//			y = f(x)
	//		��ʾ
	//

	CMatrix cMatrixOutputLayerOutput(nOutputLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());


	matrixOutputLayerOutput.CopyMatrix(cMatrixOutputLayerOutput);

}

/////////////////////////////////////////////////////////////////////////////
//	Back propagation		---->		ǰ�����(Only for Simulating)	   //
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

		// �õ�����������ֵ
		matrixDemoDataInput.CopySubMatrix (cMatrixInputLayerValue, (unsigned int)0, (unsigned int)0);

		CMatrix cMatrixTInputLayerValue = cMatrixInputLayerValue.Transpose ();

		matrixInputLayerValue.CopyMatrix (cMatrixTInputLayerValue);


	//}

	/////////////////////////////////////////////////////////////////////////
	// �õ�����������������ľ�����
	// �������:
	//		1. ��������Ŀ��Ϊ��������;
	//		2. �������������������Ŀ��Ϊ���������;
	//		3. ����Ԫ���е�ֵ��Ϊ��Ӧ��������������ľ�����:
	//		   �� 
	//				cMatrixInputLayerValue * cMatrixInputToHideWeightValue 
	//				+ cMatrixHideLayerValveValue
	//			�õ�.
	//
		
	CMatrix cMatrixExHideLayerValveValue;
	cMatrixExHideLayerValveValue.nncpyi (matrixHideLayerValveValue, matrixDemoDataInput.GetMatrixRowNumber ());

	CMatrix cMatrixHideLayerPureInput(nHideLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());

	cMatrixHideLayerPureInput = matrixInputToHideWeightValue * matrixInputLayerValue;

	cMatrixHideLayerPureInput += cMatrixExHideLayerValveValue;


	/////////////////////////////////////////////////////////////////////
	//	�����������������������
	//	�������:
	//		1. ����������y�������������x�Ĺ�ϵ���ú�����ʾ
	//			y = f(x)
	//		2. �����ά����������ľ���������ά�����
	//

	CMatrix cMatrixHideLayerOutput(nHideLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());

		cMatrixHideLayerOutput = cMatrixHideLayerPureInput.Tansig();
	
	
	matrixHideLayerOutput.CopyMatrix(cMatrixHideLayerOutput);

	/////////////////////////////////////////////////////////////////////
	//	�õ��������������ľ�����
	//	�������;
	//		1. ��������Ŀ��Ϊ���������;
	//		2. ������������������Ŀ��Ϊ���������;
	//		3. ������Ԫ�ص�ֵ��Ϊ��Ӧ�����������ľ�����:
	//			��
	//				cMatrixHideLayerOutput * cMatrixHideToOutputWeightValue
	//				+ cMatrixOutputLayerValveValue
	//			�õ�
	//

	CMatrix cMatrixExOutputLayerValveValue;
	cMatrixExOutputLayerValveValue.nncpyi (matrixOutputLayerValveValue, matrixDemoDataInput.GetMatrixRowNumber ());

	CMatrix cMatrixOutputLayerPureInput(nOutputLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());

	cMatrixOutputLayerPureInput = matrixHideToOutputWeightValue * cMatrixHideLayerOutput;
	cMatrixOutputLayerPureInput  += cMatrixExOutputLayerValveValue;


	/////////////////////////////////////////////////////////////////////
	//	����������������������
	//	�������:
	//		1. �����ά����õ������������������ľ�������ɵľ���һ��;
	//		2. ���������y��������������ù�ϵʽ
	//			y = f(x)
	//		��ʾ
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
	// ����BP�㷨����nStep�ĳ�ʼֵ
	nStep = 0.1;

	// ǰ�����
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
	//	�������������������delta����
	//	�������:
	//		1. ��������ĿΪ���������;
	//		2. ������������ĿΪ���������;
	//		3. �����е�Ԫ�ص�ֵyΪ:
	//				y = (ǰ��������������ֵ - ������������ֵ) .* f'(net)
	//
	CMatrix cMatrixTDemoOutput(matrixDemoDataInput.GetMatrixRowNumber (), nOutputLayerNumber);
		
	// �õ�����������������
	matrixDemoDataInput.CopySubMatrix (cMatrixTDemoOutput, (unsigned int)nInputLayerNumber, (unsigned int)0);

	CMatrix cMatrixDemoOutput = cMatrixTDemoOutput.Transpose ();

	// �õ����������������
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

		// ��������deltaֵ
		// ע��: �˴�'/' �� '���'!!!
		CMatrix	cMatrixOutputLayerDelta (nOutputLayerNumber, matrixDemoDataInput.GetMatrixRowNumber());
		cMatrixOutputLayerDelta = (matrixOutputLayerOutput - matrixOutputLayerOutput / matrixOutputLayerOutput) / cMatrixOutputLayerError;
	
		CMatrix cMatrixTHideToOutputWeightValue (matrixHideToOutputWeightValue.GetMatrixColNumber(), matrixHideToOutputWeightValue.GetMatrixRowNumber());
		cMatrixTHideToOutputWeightValue = matrixHideToOutputWeightValue.Transpose();

		// ���������deltaֵ
		// ע��: �˴�'/' �� '���'!!!
		CMatrix cMatrixHideLayerDelta;
		cMatrixHideLayerDelta.CopyMatrix ( (matrixHideLayerOutput - (matrixHideLayerOutput / matrixHideLayerOutput)) / ( cMatrixTHideToOutputWeightValue * cMatrixOutputLayerDelta) );
				
		// �����µ�����㵽�������Ȩֵ
		CMatrix cMatrixNewInputToHideWeight (matrixInputToHideWeightValue.GetMatrixRowNumber (), matrixInputToHideWeightValue.GetMatrixColNumber ());
		// ������µ�������ķ�ֵ
		CMatrix cMatrixNewHideLayerValve (nHideLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());
		// �����µ������㵽������Ȩֵ
		CMatrix cMatrixNewHideToOutputWeight (matrixHideToOutputWeightValue.GetMatrixRowNumber (), matrixHideToOutputWeightValue.GetMatrixColNumber ());
		// �����µ������ķ�ֵ
		CMatrix cMatrixNewOutputLayerValve (nOutputLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());
		// �����µ�������
		CMatrix cMatrixNewOutputLayerError(nOutputLayerNumber, matrixDemoDataInput.GetMatrixRowNumber ());


		// Ȩֵ�ͷ�ֵ����
		cMatrixNewHideToOutputWeight = cMatrixOutputLayerDelta * (matrixHideLayerOutput.Transpose ()) * (nStep);
		cMatrixNewOutputLayerValve = cMatrixOutputLayerDelta;

		cMatrixNewInputToHideWeight = cMatrixHideLayerDelta * (matrixInputLayerValue.Transpose ()) * (nStep);
		cMatrixNewHideLayerValve = cMatrixHideLayerDelta;

		// ��ֵ
		matrixInputToHideWeightValue += cMatrixNewInputToHideWeight;

		CMatrix cMatrixExHideLayerValveValue;
		cMatrixExHideLayerValveValue.nncpyi (matrixHideLayerValveValue, matrixDemoDataInput.GetMatrixRowNumber ());
		cMatrixExHideLayerValveValue += cMatrixNewHideLayerValve;

		matrixHideToOutputWeightValue += cMatrixNewHideToOutputWeight;

		CMatrix cMatrixExOutputLayerValveValue;
		cMatrixExOutputLayerValveValue.nncpyi (matrixOutputLayerValveValue, matrixDemoDataInput.GetMatrixRowNumber ());
		cMatrixExOutputLayerValveValue += cMatrixNewOutputLayerValve;

		// ǰ�����
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

		// ��ʾ���ݺͳ�������״̬
		nSystemError = nSystemErrorOld;
		nTrainTimes = nLoopTimes;


		// ��ʾϵͳ���
	/*	CString	strSystemError;
		strSystemError.Format ("%lf", nSystemError);
		LPCTSTR	lpstrSystemError = (LPCTSTR)strSystemError;
		HWND	hwnd = ::GetDlgItem (hWnd, ID_SYSTEM_ERROR);
		::SetWindowText (hwnd, lpstrSystemError);*/
	
		// ��ʾѵ������
		/*CString strTrainTimes;
		strTrainTimes.Format ("%u", nTrainTimes + 1);
		LPCTSTR lpstrTrainTimes = (LPCTSTR)strTrainTimes;
		hwnd = ::GetDlgItem (hWnd, ID_TRAIN_TIMES);
		::SetWindowText (hwnd, lpstrTrainTimes);*/

	}// End the "for" loop
	
	return true;
}



