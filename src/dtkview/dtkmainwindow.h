#ifndef DTKMAINWINDOW_H
#define DTKMAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class DTKMainWindow;
}

class DTKMainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit DTKMainWindow(QWidget *parent = 0);
    ~DTKMainWindow();

public slots:
    void updateStatusBar(void);

private slots:
    void on_resetPushButton_clicked();

private:
    Ui::DTKMainWindow *ui;
};

#endif // DTKMAINWINDOW_H
