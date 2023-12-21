// widgetstream.cpp

#include "widgetstream.h"
#include <QApplication>
WidgetStream::WidgetStream(std::ostream &stream, QTextEdit *text_edit)
    : m_stream(stream), m_text_edit(text_edit) {
    m_old_buf = stream.rdbuf();
    stream.rdbuf(this);
}

WidgetStream::~WidgetStream() {
    m_stream.rdbuf(m_old_buf);
}

std::streambuf::int_type WidgetStream::overflow(int_type v) {
    if (v == '\n') {
        log_window_update();
    }
    return v;
}
int WidgetStream::sync() {
    log_window_update(m_string);
    m_string.clear();
    return 0; // Returning zero indicates success
}
std::streamsize WidgetStream::xsputn(const char *p, std::streamsize n) {
    m_string.append(p, p + n);

    size_t pos = 0;
    while ((pos = m_string.find('\n')) != std::string::npos) {
        std::string tmp(m_string.begin(), m_string.begin() + pos);
        log_window_update(tmp);
        m_string.erase(m_string.begin(), m_string.begin() + pos + 1);
    }

    return n;
}
void WidgetStream::log_window_update(const std::string &s) {
    if (!s.empty()) {
        m_text_edit->append(QString::fromStdString(s));
        QApplication::processEvents(); // Process the event loop to update the GUI
    }
}


