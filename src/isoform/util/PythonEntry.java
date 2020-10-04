package isoform.util;

import py4j.GatewayServer;

public class PythonEntry {

    public static void main(String[] args) {
    	int port=Integer.parseInt(args[0]);
    	System.out.println("Java: starting server on "+port);
        GatewayServer gatewayServer = new GatewayServer(null, port);
        gatewayServer.start();
        System.out.println("Gateway Server Started");
    }

}