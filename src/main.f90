PROGRAM run_me
    IMPLICIT NONE

    ! 1. Skriv til terminalvinduet
    PRINT *, "Hello, World! Modern Fortran er i gang på PC-en!"

    ! 2. Skriv til en tekstfil (nyttig for å bekrefte filoperasjoner)
    OPEN(UNIT=10, FILE="output.txt", STATUS="REPLACE", ACTION="WRITE")
    WRITE(10, '(A)') "Dette er skrevet til fil fra en Fortran-rutine."
    CLOSE(10)

    PRINT *, "Tekstfilen 'output.txt' ble opprettet og skrevet til."
END PROGRAM run_me