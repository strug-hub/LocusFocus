# Development image
FROM python:3.7
WORKDIR /usr/src/app
COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt
COPY . .
EXPOSE 5000
ENV FLASK_ENV development
ENV TZ America/Toronto
ENTRYPOINT [ "python3", "-m" , "flask", "run" ]
CMD [ "--host=0.0.0.0"]
