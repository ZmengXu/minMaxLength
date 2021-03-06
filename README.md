# minMaxLength
This is a function object which can return the maximum and minimum length of a field reaching a criterion value. For example, you can use it to get the temporal evolution of the lift-off length, e.g., setting the fieldName as OH, criterion as a mass fraction of 0.001. You can also use it to get the vapor penetration length, e.g., setting the fieldName as Z, criterion as a Z = 0.001.


```
    minMaxLength1
    {
        type minMaxLength;
        libs ("libSJfieldFunctionObjects.so");
        fields
        (
            Z
            {
                criterion       0.001;
                threshold       0.0001;
            }
            OH
            {
                criterion       0.0001;//2% of 0.005
            }
        );
        position        ( 0 0 0 );
        direction       ( 0 0 1 );
    }
```
