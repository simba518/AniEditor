#include <stdlib.h>
#include <QtGui/QApplication>
#include <Log.h>
#include <QInputEventRecorderCmd.h>
#include "AniEditMainWin.h"
using namespace UTILITY;

#ifdef WIN32
const string TEST_UI_INIF= "D:\\lsw\\AnimationEditor\\TestCase\\TestData\\ani_ui_test_windows.ini";
#else
const string TEST_UI_INIF= "/home/simba/Workspace/AnimationEditor/TestCase/TestData/ani_ui_test.ini";
#endif

int main(int argc, char *argv[]){

  QApplication application(argc,argv);
  QMainWindow w;
  ANI_EDIT_UI::AniEditMainWin m_mainw("",&w);

  // auto record and replay
  pQInputEventRecorderCmd recorder;
  if(argc >= 2){
	recorder = pQInputEventRecorderCmd(new QInputEventRecorderCmd(&w));
	if (2 == argc)
	  recorder->setCmd(argv[argc-1]);
	else if(3 == argc)
	  recorder->setCmd(argv[argc-2],argv[argc-1]);
  }

  if(!m_mainw.loadInitFile(TEST_UI_INIF))
	ERROR_LOG("failed to load the init file: "<<TEST_UI_INIF);

  m_mainw.show();
  return  application.exec();
}
