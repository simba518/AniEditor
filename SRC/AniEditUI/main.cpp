#include <stdlib.h>
#include <QtGui/QApplication>
#include <Log.h>
#include "AniEditMainWin.h"
using namespace UTILITY;

int main(int argc, char *argv[]){

  const string home_dir = string(getenv("HOME"));
  QApplication application(argc,argv);
  const string main_inf = home_dir + "/Workspace/PathControl/MainWindow.ini";
  LSW_ANI_EDIT_UI::AniEditMainWin m_mainw(main_inf);
  m_mainw.show();
  return  application.exec();
}
