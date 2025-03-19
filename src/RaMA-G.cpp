// RaMA-G.cpp: 定义应用程序的入口点。
//

#include <iostream>

#include "config.hpp"



int main(int argc, char** argv)
{
	spdlog::init_thread_pool(8192, 1);
	setupLogger();

	CLI::App app{ "RaMA-G: A High-performance Genome Alignment Tool" };



	CommonArgs common_args;

	addCommonOptions(&app, common_args);


	CLI11_PARSE(app, argc, argv);

	spdlog::info("Hello, World!");

	spdlog::info("Hello, World!");

}