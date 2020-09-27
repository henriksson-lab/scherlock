from py4j.java_gateway import JavaGateway
gateway = JavaGateway()                        # connect to the JVM
java_object = gateway.jvm.mypackage.MyClass()  # invoke constructor
other_object = java_object.doThat()
other_object.doThis(1,'abc')
gateway.jvm.java.lang.System.out.println('Hello World!') # call a static method


