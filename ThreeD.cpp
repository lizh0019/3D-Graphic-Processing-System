 // ThreeD.cpp : implementation file
//

#include "stdafx.h"
#include "kettle.h"
#include "ThreeD.h"
// Standard include files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <time.h>
#include <MATH.H>
#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
/////////////////////////////////////////////////////////////////////////////
// ThreeD

IMPLEMENT_DYNCREATE(ThreeD, CDocument)

ThreeD::ThreeD()
{
	opensuccess=FALSE;
	max=0;
	center[0]=center[1]=center[2]=0;
	m_FileName="";
	pointnumber=0;
	int i,j;
	for(i=0;i<MAXPOINTNUMBER;i++){
		for(j=0;j<3;j++){
			data[i][j]=0;
			data1[i][j]=0;
		}
	}
}
ThreeD::ThreeD(ThreeD &thd){
	int i,j;
	for(i=0;i<MAXPOINTNUMBER;i++){
		for(j=0;j<3;j++){
			data[i][j]=thd.data[i][j];
			data1[i][j]=thd.data1[i][j];
		}
	}
	m_FileName=thd.m_FileName;
	opensuccess=thd.opensuccess;
	pointnumber=thd.pointnumber;
	max=thd.max;
	for(i=0;i<3;i++){
		center[i]=thd.center[i];
	}
}
BOOL ThreeD::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;
	return TRUE;
}

ThreeD::~ThreeD()
{
}


BEGIN_MESSAGE_MAP(ThreeD, CDocument)
	//{{AFX_MSG_MAP(ThreeD)
		// NOTE - the ClassWizard will add and remove mapping macros here.
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// ThreeD diagnostics

#ifdef _DEBUG
void ThreeD::AssertValid() const
{
	CDocument::AssertValid();
}

void ThreeD::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// ThreeD serialization

void ThreeD::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: add storing code here
	}
	else
	{
		// TODO: add loading code here
	}
}

/////////////////////////////////////////////////////////////////////////////
// ThreeD commands

BOOL ThreeD::Read3D_RMA()
{
	CFile fin;
	CFileException e;
	short x,y,z;
	short cnt=0;
	int stripe=0;
	int numperline;
	long i;int j; 
    CString  showxyz;
    int size=0;
	//pointnumber=0;
	CFileDialog openxyz(TRUE,NULL,NULL,OFN_HIDEREADONLY,"XYZ File(*.xyz)|*.xyz|");
	openxyz.m_ofn.lpstrTitle="打开xyz文件";
	if(!(openxyz.DoModal()==IDOK)){
         AfxMessageBox("打开xyz文件错误！",MB_ICONSTOP|MB_OK);
		 return FALSE;
	}
	m_FileName=openxyz.GetFileName();
	if (!fin.Open(m_FileName,CFile::modeRead,&e)) return FALSE;
   /* Read data from file */
   while (fin.Read(&cnt,sizeof(short)) == sizeof(short)){	// *********** Read raw data *************
		/* New stripe */
		numperline=0;
		for (i=0; i < cnt; i++) {
			if (fin.Read(&x, sizeof(short)) <1) return FALSE;
			if (fin.Read(&y, sizeof(short)) <1) return FALSE;
			if (fin.Read(&z, sizeof(short)) <1) return FALSE;
		    //TRACE("(%6d,%6d,%6d)\n",x,y,z);
            data[pointnumber][0]=x;
            data[pointnumber][1]=y;
            data[pointnumber][2]=z;
			pointnumber++;
		    numperline++;
		}
		//TRACE("stripe%4d: %6d\n",stripe,numperline);
		stripe++;
   }
   //TRACE("*** stripes ***:  %6d\n",stripe);
	for(j=0;j<3;j++){/*find the the center*/
		for(i=0;i<pointnumber;i++){
			center[j]+=data[i][j];
		}
		center[j]/=pointnumber;
	}
	for(i=0;i<pointnumber;i++){/*translation and find the max element*/
		for(j=0;j<3;j++){
			data1[i][j]=data[i][j]-center[j];
			if(max<fabs(data1[i][j]))
				max=fabs(data1[i][j]);
		}
	}
	for(i=0;i<pointnumber;i++){/*normalization and store data in data1[MAXPOINTNUMBER][3]*/
		for(j=0;j<3;j++){
			data1[i][j]/=(float)max;
		}
	}
   fin.Close();
   opensuccess=true;
   return pointnumber;
}

BOOL ThreeD::Save3D_RMAData(int pointnumber)
{
	if(!opensuccess){
		AfxMessageBox("3D model has not yet be loaded");
		return FALSE;
	}
	CFileException e;
	char BASED_CODE szFilter[]="3D Model Vertices Data(*.txt)|*.txt||";
	CFileDialog m_ldFile(FALSE,"*.txt",m_FileName,OFN_HIDEREADONLY|OFN_OVERWRITEPROMPT,szFilter);
	m_ldFile.m_ofn.lpstrTitle="存储xyz文件";
	if(m_ldFile.DoModal()!=IDOK)
		return FALSE;
	m_FileName=m_ldFile.GetPathName();
	//write file
	FILE *fpout=NULL;
	CString dataout="",temp="";
//	if (!fpout.Open(m_FileName,CFile::modeWrite,&e)) return FALSE;
	if (!(fpout=fopen(m_FileName,"w"))) {
		temp.Format("Unable to open file %s",m_FileName);
		AfxMessageBox(temp);
		return FALSE;
	}
	//write pixel information
	for (long i=0;i<pointnumber;i++){
		temp.Format("%d %d %d\n",data[i][0],data[i][1],data[i][2]);
		dataout+=temp;
	}
	//fpout->Write(dataout,pointnumber*3);
	//fpout->Close();
	fputs(dataout,fpout);
	fclose(fpout);
	return TRUE;
}


