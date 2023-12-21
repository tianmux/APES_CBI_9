// widgetstream.h

#ifndef WIDGETSTREAM_H
#define WIDGETSTREAM_H

#include <iostream>
#include <streambuf>
#include <string>
#include <QTextEdit>

class WidgetStream : public std::streambuf {
public:
    WidgetStream(std::ostream &stream, QTextEdit *text_edit);
    ~WidgetStream();

protected:
    virtual int_type overflow(int_type v);
    virtual std::streamsize xsputn(const char *p, std::streamsize n);
    int sync();
private:
    void log_window_update(const std::string &s = "");

    std::ostream &m_stream;
    std::streambuf *m_old_buf;
    QTextEdit *m_text_edit;
    std::string m_string;
};

#endif // WIDGETSTREAM_H
