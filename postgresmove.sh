#!/bin/sh

for i in */*.py
do
ed $i <<END
/import db
d
.
/import reportlib
a

try:
    if os.environ['DB_TYPE'] == 'postgres':
        import pg_db
        db = pg_db
        db.setTrace()
        db.setAutoTranslateBE()
    else:
        import db
except:
    import db

.
w
q
END
done

