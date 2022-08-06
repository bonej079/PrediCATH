# Import smtplib for the actual sending function
import smtplib


class EmailManager:
    @staticmethod
    def send_message(to: str, title: str, message: str):
        try:
            from email.mime.text import MIMEText
            import socket
            message += f"\n\nMessage sent from {socket.gethostname()}"
            msg = MIMEText(message)

            # me == the sender's email address
            # you == the recipient's email address
            msg['Subject'] = title
            msg['From'] = 'joseph.bonello@um.edu.mt'
            msg['To'] = to

            # Send the message via our own SMTP server.
            s = smtplib.SMTP('smtp.gmail.com', 587, timeout=30)
            s.ehlo()
            s.starttls()
            s.login('jbphdcode@gmail.com', 'efvztygnxrwambsk')  # '8Kh%7eonhBlz')
            s.send_message(msg)
            s.quit()
        except:
            print("Sending the message failed!")


def main():
    email = EmailManager()
    email.send_message('joseph.bonello@gmail.com', 'Python Utility Test', 'Try succeeded')
    print("Finished")


if __name__ == '__main__':
    main()
