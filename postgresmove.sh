#!/bin/sh

for i in postgres/*.py
do
ed $i <<END
/import db
d
.
/import reportlib
a

if os.environ['DB_TYPE'] == 'postgres':
    import pg_db
    db = pg_db
    db.setTrace()
    db.setAutoTranslateBE()
else:
    import db

.
w
q
END
done

