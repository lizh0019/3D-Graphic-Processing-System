; CLW file contains information for the MFC ClassWizard

[General Info]
Version=1
LastClass=CKettleDoc
LastTemplate=CDocument
NewFileInclude1=#include "stdafx.h"
NewFileInclude2=#include "kettle.h"
LastPage=0

ClassCount=7
Class1=CKettleApp
Class2=CKettleDoc
Class3=CKettleView
Class4=CMainFrame

ResourceCount=4
Resource1=IDR_KETTLETYPE
Resource2=IDD_ABOUTBOX
Class5=CChildFrame
Class6=CAboutDlg
Class7=ThreeD
Resource3=IDR_MAINFRAME
Resource4=IDD_DIALOG1

[CLS:CKettleApp]
Type=0
HeaderFile=kettle.h
ImplementationFile=kettle.cpp
Filter=N

[CLS:CKettleDoc]
Type=0
HeaderFile=kettleDoc.h
ImplementationFile=kettleDoc.cpp
Filter=N
BaseClass=CDocument
VirtualFilter=DC
LastObject=ID_FILE_OPEN

[CLS:CKettleView]
Type=0
HeaderFile=kettleView.h
ImplementationFile=kettleView.cpp
Filter=C
LastObject=IDM_PROCESS
BaseClass=CView
VirtualFilter=VWC


[CLS:CMainFrame]
Type=0
HeaderFile=MainFrm.h
ImplementationFile=MainFrm.cpp
Filter=T
BaseClass=CMDIFrameWnd
VirtualFilter=fWC


[CLS:CChildFrame]
Type=0
HeaderFile=ChildFrm.h
ImplementationFile=ChildFrm.cpp
Filter=M


[CLS:CAboutDlg]
Type=0
HeaderFile=kettle.cpp
ImplementationFile=kettle.cpp
Filter=D

[DLG:IDD_ABOUTBOX]
Type=1
Class=CAboutDlg
ControlCount=4
Control1=IDC_STATIC,static,1342177283
Control2=IDC_STATIC,static,1342308480
Control3=IDC_STATIC,static,1342308352
Control4=IDOK,button,1342373889

[MNU:IDR_MAINFRAME]
Type=1
Class=CMainFrame
Command1=ID_FILE_NEW
Command2=ID_FILE_OPEN
Command3=ID_FILE_PRINT_SETUP
Command4=ID_FILE_MRU_FILE1
Command5=ID_APP_EXIT
Command6=ID_VIEW_TOOLBAR
Command7=ID_VIEW_STATUS_BAR
Command8=ID_APP_ABOUT
CommandCount=8

[TB:IDR_MAINFRAME]
Type=1
Class=CMainFrame
Command1=ID_FILE_NEW
Command2=ID_FILE_OPEN
Command3=ID_FILE_SAVE
Command4=ID_EDIT_CUT
Command5=ID_EDIT_COPY
Command6=ID_EDIT_PASTE
Command7=ID_FILE_PRINT
Command8=ID_APP_ABOUT
CommandCount=8

[MNU:IDR_KETTLETYPE]
Type=1
Class=CKettleView
Command1=ID_FILE_NEW
Command2=ID_FILE_OPEN
Command3=ID_FILE_CLOSE
Command4=ID_FILE_SAVE
Command5=ID_FILE_SAVE_AS
Command6=ID_FILE_PRINT
Command7=ID_FILE_PRINT_PREVIEW
Command8=ID_FILE_PRINT_SETUP
Command9=ID_FILE_MRU_FILE1
Command10=ID_APP_EXIT
Command11=ID_OPEN_3DRMA
Command12=ID_EDIT_UNDO
Command13=ID_EDIT_CUT
Command14=ID_EDIT_COPY
Command15=ID_EDIT_PASTE
Command16=ID_VIEW_TOOLBAR
Command17=ID_VIEW_STATUS_BAR
Command18=ID_WINDOW_NEW
Command19=ID_WINDOW_CASCADE
Command20=ID_WINDOW_TILE_HORZ
Command21=ID_WINDOW_ARRANGE
Command22=IDM_LINES
Command23=IDM_POINTS
Command24=IDM_SHADE
Command25=IDM_TEXTURED
Command26=IDM_ROTATE
Command27=IDM_TRIANGLE
Command28=IDM_PROCESS
Command29=ID_APP_ABOUT
CommandCount=29

[ACL:IDR_MAINFRAME]
Type=1
Class=CMainFrame
Command1=ID_FILE_NEW
Command2=ID_FILE_OPEN
Command3=ID_FILE_SAVE
Command4=ID_FILE_PRINT
Command5=ID_EDIT_UNDO
Command6=ID_EDIT_CUT
Command7=ID_EDIT_COPY
Command8=ID_EDIT_PASTE
Command9=ID_EDIT_UNDO
Command10=ID_EDIT_CUT
Command11=ID_EDIT_COPY
Command12=ID_EDIT_PASTE
Command13=ID_NEXT_PANE
Command14=ID_PREV_PANE
CommandCount=14

[CLS:ThreeD]
Type=0
HeaderFile=ThreeD.h
ImplementationFile=ThreeD.cpp
BaseClass=CDocument
Filter=N

[DLG:IDD_DIALOG1]
Type=1
Class=?
ControlCount=2
Control1=IDOK,button,1342242817
Control2=IDCANCEL,button,1342242816

